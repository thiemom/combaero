"""A minimal in-process HTTP client for the GUI's FastAPI app.

`fastapi.testclient` needs `httpx`, which is not a dependency here. This
drives the ASGI app directly with one request, so a test goes through the
whole HTTP layer the browser does: routing, JSON body parsing, request
validation (422), the route body, `response_model` serialization and the
JSON encoder. Calling a route coroutine directly skips all of that except
the route body.
"""

from __future__ import annotations

import asyncio
import json
import time
from dataclasses import dataclass, field
from typing import Any


@dataclass
class Response:
    status_code: int
    headers: dict[str, str] = field(default_factory=dict)
    body: bytes = b""
    #: Seconds from the request to the response start, as the browser sees it
    #: (asyncio.run may wait longer, for worker threads, before returning).
    elapsed: float = 0.0

    @property
    def text(self) -> str:
        return self.body.decode("utf-8", errors="replace")

    def json(self) -> Any:
        return json.loads(self.body)


async def _call(app: Any, method: str, path: str, body: bytes, headers: list) -> Response:
    scope = {
        "type": "http",
        "asgi": {"version": "3.0"},
        "http_version": "1.1",
        "method": method,
        "scheme": "http",
        "path": path,
        "raw_path": path.encode(),
        "query_string": b"",
        "root_path": "",
        "headers": headers,
        "client": ("127.0.0.1", 50000),
        "server": ("127.0.0.1", 8000),
    }
    sent = False

    async def receive() -> dict:
        nonlocal sent
        if not sent:
            sent = True
            return {"type": "http.request", "body": body, "more_body": False}
        # The app only reads past the body to watch for a disconnect; never
        # report one, so a slow route is not cancelled by the client side.
        await asyncio.Event().wait()
        return {"type": "http.disconnect"}

    response = Response(status_code=0)
    t0 = time.perf_counter()
    chunks: list[bytes] = []

    async def send(message: dict) -> None:
        if message["type"] == "http.response.start":
            response.status_code = message["status"]
            response.elapsed = time.perf_counter() - t0
            response.headers = {k.decode().lower(): v.decode() for k, v in message["headers"]}
        elif message["type"] == "http.response.body":
            chunks.append(message.get("body", b""))

    try:
        await app(scope, receive, send)
    except Exception:
        # Starlette's ServerErrorMiddleware sends its plain 500 and then
        # re-raises for the server to log; the response is what the browser
        # sees, so keep that and let a test assert on it.
        if response.status_code == 0:
            raise
    response.body = b"".join(chunks)
    return response


def request(
    app: Any,
    method: str,
    path: str,
    json_body: Any = None,
    raw_body: bytes | None = None,
    host: str = "127.0.0.1:8000",
) -> Response:
    """Send one request to `app` and return the complete response."""
    if raw_body is None:
        raw_body = b"" if json_body is None else json.dumps(json_body).encode()
    headers = [(b"host", host.encode())]
    if raw_body:
        headers.append((b"content-type", b"application/json"))
        headers.append((b"content-length", str(len(raw_body)).encode()))
    return asyncio.run(_call(app, method, path, raw_body, headers))
