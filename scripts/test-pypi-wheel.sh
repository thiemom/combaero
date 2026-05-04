#!/usr/bin/env bash
# Verify the combaero-gui wheel end-to-end before publishing to PyPI.
#
# Mirrors what the CI publish-gui.yml verify job does, but runs locally.
# All build artefacts and the test venv are written to .package_tests/ (gitignored).
# The combaero core wheel (C++ extension) is built to dist/ as normal.
#
# Usage:
#   scripts/test-pypi-wheel.sh [--skip-frontend] [--skip-combaero-build]
#
# Options:
#   --skip-frontend        Skip pnpm build if gui/frontend/dist/ already exists.
#   --skip-combaero-build  Skip building the combaero core wheel (uv build --package combaero).
#                          If dist/combaero-*.whl already exists it is used for the test;
#                          otherwise combaero is pulled from PyPI. Use this flag for fast
#                          iteration when only the GUI changed and combaero is unchanged.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_DIR="$ROOT_DIR/.package_tests"
DIST_DIR="$TEST_DIR/dist"
VENV_DIR="$TEST_DIR/venv"
PORT=8000

cd "$ROOT_DIR"

# ── Colours ──────────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m'

step()  { echo -e "\n${BLUE}${BOLD}── $* ──${NC}"; }
ok()    { echo -e "  ${GREEN}OK${NC}  $*"; }
info()  { echo -e "  ${YELLOW}...${NC} $*"; }
fail()  { echo -e "  ${RED}FAIL${NC} $*" >&2; }

# ── CLI flags ─────────────────────────────────────────────────────────────────
SKIP_FRONTEND=false
SKIP_COMBAERO_BUILD=false
for arg in "$@"; do
    case "$arg" in
        --skip-frontend)       SKIP_FRONTEND=true ;;
        --skip-combaero-build) SKIP_COMBAERO_BUILD=true ;;
        --help|-h)
            echo "Usage: $0 [--skip-frontend] [--skip-combaero-build]"
            echo "  --skip-frontend        Reuse existing gui/frontend/dist/ instead of rebuilding."
            echo "  --skip-combaero-build  Skip uv build --package combaero (uses dist/ wheel or PyPI)."
            exit 0
            ;;
    esac
done

# ── Cleanup trap ──────────────────────────────────────────────────────────────
SERVER_PID=""
cleanup() {
    if [ -n "$SERVER_PID" ]; then
        kill "$SERVER_PID" 2>/dev/null || true
        SERVER_PID=""
    fi
}
trap cleanup EXIT INT TERM

echo ""
echo -e "${BOLD}combaero-gui wheel verification${NC}"
echo "  Test artefacts → $TEST_DIR"
echo "  Port           → $PORT"

# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Build frontend
# ─────────────────────────────────────────────────────────────────────────────
step "Step 1  Build frontend"

if [ "$SKIP_FRONTEND" = true ] && [ -d "gui/frontend/dist" ]; then
    info "Skipping pnpm build (--skip-frontend, dist/ exists)"
else
    info "pnpm install + build..."
    (cd gui/frontend && pnpm install --frozen-lockfile && pnpm run build)
    ok "Frontend built → gui/frontend/dist/"
fi

if [ ! -f "gui/frontend/dist/index.html" ]; then
    fail "gui/frontend/dist/index.html not found — frontend build may have failed"
    exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Build combaero core wheel
# ─────────────────────────────────────────────────────────────────────────────
step "Step 2  Build combaero core wheel"

if [ "$SKIP_COMBAERO_BUILD" = true ]; then
    COMBAERO_WHEEL=$(ls "$ROOT_DIR/dist"/combaero-*.whl 2>/dev/null | head -1 || true)
    if [ -n "$COMBAERO_WHEEL" ]; then
        info "Skipping build (--skip-combaero-build), using: $(basename "$COMBAERO_WHEEL")"
    else
        info "Skipping build (--skip-combaero-build), no local wheel found — will pull from PyPI"
    fi
else
    info "uv build --package combaero ..."
    uv build --package combaero
    COMBAERO_WHEEL=$(ls "$ROOT_DIR/dist"/combaero-*.whl 2>/dev/null | head -1 || true)
    if [ -z "$COMBAERO_WHEEL" ]; then
        fail "No combaero wheel found in dist/ after build"
        exit 1
    fi
    ok "Built: $(basename "$COMBAERO_WHEEL")"
fi

# ─────────────────────────────────────────────────────────────────────────────
# Step 3: Build combaero-gui wheel
# ─────────────────────────────────────────────────────────────────────────────
step "Step 3  Build combaero-gui wheel"

rm -rf "$DIST_DIR"
info "uv build --package combaero-gui ..."
uv build --package combaero-gui --out-dir "$DIST_DIR"

GUI_WHEEL=$(ls "$DIST_DIR"/combaero_gui-*.whl 2>/dev/null | head -1 || true)
if [ -z "$GUI_WHEEL" ]; then
    fail "No .whl found in $DIST_DIR"
    exit 1
fi
ok "Built: $(basename "$GUI_WHEEL")"

# ─────────────────────────────────────────────────────────────────────────────
# Step 4: Inspect wheel contents
# ─────────────────────────────────────────────────────────────────────────────
step "Step 4  Inspect wheel contents"

python3 - "$GUI_WHEEL" <<'PYEOF'
import zipfile, sys

wheel = sys.argv[1]
with zipfile.ZipFile(wheel) as z:
    names = z.namelist()

html = [n for n in names if n.endswith("index.html")]
js   = [n for n in names if n.endswith(".js")  and "/assets/" in n]
css  = [n for n in names if n.endswith(".css") and "/assets/" in n]

def _short(paths, limit=3):
    shown = [p.split("/", 1)[-1] for p in paths[:limit]]
    rest  = len(paths) - limit
    suffix = f" … (+{rest} more)" if rest > 0 else ""
    return ", ".join(shown) + suffix

print(f"  index.html : {_short(html)}")
print(f"  JS  assets : {_short(js)}")
print(f"  CSS assets : {_short(css)}")
print(f"  total      : {len(names)} files in wheel")

errors = []
if not html:
    errors.append("index.html is missing from wheel")
if not js:
    errors.append(
        "No JS files found under assets/ — frontend/dist/assets/ was not bundled.\n"
        "  This is the white-page MIME bug: the /assets StaticFiles mount will be\n"
        "  skipped at startup and the SPA catch-all will serve HTML for every asset."
    )

if errors:
    for e in errors:
        print(f"\nFAIL: {e}", file=sys.stderr)
    sys.exit(1)
PYEOF

ok "Wheel contents verified"

# ─────────────────────────────────────────────────────────────────────────────
# Step 5: Install in a clean venv and smoke-test the server
# ─────────────────────────────────────────────────────────────────────────────
step "Step 5  Install wheel in a clean venv and smoke-test"

if lsof -ti:"$PORT" >/dev/null 2>&1; then
    fail "Port $PORT is already in use — stop the process and retry."
    exit 1
fi

info "Creating isolated venv in $VENV_DIR ..."
rm -rf "$VENV_DIR"
uv venv "$VENV_DIR" --python 3.12 --quiet

# Resolve COMBAERO_WHEEL in case --skip-combaero-build set it (or left it empty)
COMBAERO_WHEEL="${COMBAERO_WHEEL:-}"
COMBAERO_WHEEL=$(ls "$ROOT_DIR/dist"/combaero-*.whl 2>/dev/null | head -1 || true)

if [ -n "$COMBAERO_WHEEL" ]; then
    info "Installing: $(basename "$COMBAERO_WHEEL") + $(basename "$GUI_WHEEL") ..."
    # Install both together so uv resolves all deps (fastapi, uvicorn, etc.) in one pass
    # and does not pull combaero from PyPI.
    uv pip install --python "$VENV_DIR" "$COMBAERO_WHEEL" "$GUI_WHEEL" --quiet
else
    info "Installing: $(basename "$GUI_WHEEL") (combaero from PyPI) ..."
    uv pip install --python "$VENV_DIR" "$GUI_WHEEL" --quiet
fi
ok "Installation complete"

info "Starting combaero-gui on port $PORT ..."
"$VENV_DIR/bin/combaero-gui" &
SERVER_PID=$!

# Wait for server to be ready — retry loop avoids a fixed sleep
READY=false
for i in $(seq 1 15); do
    if curl -fsS "http://localhost:$PORT/health" >/dev/null 2>&1; then
        READY=true
        break
    fi
    sleep 1
done

if [ "$READY" = false ]; then
    fail "Server did not become ready within 15 s — check for import errors above"
    exit 1
fi
ok "Server is up"

# ── /health ───────────────────────────────────────────────────────────────────
HEALTH=$(curl -fsS "http://localhost:$PORT/health")
echo "$HEALTH" | grep -q '"ok"' || { fail "/health returned unexpected: $HEALTH"; exit 1; }
ok "/health → $HEALTH"

# ── / must serve HTML, not the JSON API fallback ──────────────────────────────
# Use GET — modern Starlette does not auto-generate HEAD for @app.get() routes.
ROOT_BODY=$(curl -s "http://localhost:$PORT/")
echo "$ROOT_BODY" | grep -qi "<!doctype html>" || {
    fail "/ did not return HTML — got JSON API fallback instead"
    fail "(This means _FRONTEND_DIST was not found at runtime.)"
    exit 1
}
ok "/ → HTML (index.html confirmed)"

# ── JS asset must NOT be served as text/html ──────────────────────────────────
JS_PATH=$(echo "$ROOT_BODY" | grep -o 'src="[^"]*\.js"' | head -1 | sed 's/src="//;s/"//' || true)

if [ -z "$JS_PATH" ]; then
    fail "Could not find a <script src=\"*.js\"> in the served index.html"
    fail "Inspect: curl http://localhost:$PORT/"
    exit 1
fi
info "JS asset path detected: $JS_PATH"

JS_CT=$(curl -s -o /dev/null -w '%{content_type}' "http://localhost:$PORT${JS_PATH}")
if echo "$JS_CT" | grep -qi "text/html"; then
    fail "JS asset served as text/html — the MIME bug is present!"
    fail "  URL          : http://localhost:$PORT${JS_PATH}"
    fail "  content-type : $JS_CT"
    fail "  This means the /assets StaticFiles mount was not registered."
    fail "  Check that frontend/dist/assets/ is included in the wheel."
    exit 1
fi
ok "JS asset content-type: $JS_CT"

kill "$SERVER_PID" 2>/dev/null || true
SERVER_PID=""

# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo -e "${GREEN}${BOLD}All checks passed — wheel is safe to publish.${NC}"
echo ""
