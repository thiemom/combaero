# Changelog fragments

One file per user-visible change. These files replace editing `CHANGELOG.md`'s `[Unreleased]` block, which every open PR wrote at the same spot and so conflicted on.

## Adding an entry

Create `changelog.d/<id>.<category>.md`:

- **`<id>`:** the issue or PR number (`496`), or a short slug when there is neither (`changelog-fragments`).
- **`<category>`:** one of `added`, `changed`, `deprecated`, `removed`, `fixed`, `security`, `documentation`. These are the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) sections.
- **Contents:** one or more bullets in the usual style. Start with `- **Headline.** ...`, nest details as sub-bullets, and use no headings, since the category comes from the file name.

```markdown
- **The wall relay's pressure column (#496).** A wall's heat moves with ...
  - **Before:** ...
  - **Now:** ...
```

Two changes under one issue take two files: `481-merge.fixed.md`, `481-clamp.fixed.md`.

## Commands

```bash
uv run scripts/changelog.py check            # validate (also a pre-commit hook)
uv run scripts/changelog.py preview          # the release section as it would be written
uv run scripts/changelog.py release 0.7.0    # fold into CHANGELOG.md and delete the fragments
```

`release` merges the fragments with anything still under `[Unreleased]` in `CHANGELOG.md`. It writes one section per category in Keep a Changelog order, newest issue first, dates the block and moves the compare links on. Then commit the result, open the release PR, and tag as `CLAUDE.md` describes.
