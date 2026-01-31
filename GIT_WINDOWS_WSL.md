# Working from Windows and WSL — recommended Git setup

This repository enforces LF line endings via `.gitattributes`. Follow the short recommendations below so Windows and WSL users interoperate cleanly.

- Rely on `.gitattributes` (already added) to normalize line endings.
- Recommended per-environment Git configuration:
  - On Windows (Git for Windows / PowerShell / CMD):
    - `git config --global core.autocrlf true`
    - `git config --global core.filemode false`
  - On WSL / Linux:
    - `git config --global core.autocrlf input`
    - `git config --global core.filemode false`
- Alternatively, to avoid automatic conversions, set `core.autocrlf=false` on all environments and rely entirely on `.gitattributes`.

Commands to apply locally (example):

```
git add .gitattributes GIT_WINDOWS_WSL.md
git commit -m "Add .gitattributes and Windows/WSL Git instructions"
git push origin main
```

If you want, I can commit these files for you and push, or just leave them uncommitted for you to review.
