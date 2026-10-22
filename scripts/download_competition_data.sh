#!/usr/bin/env bash
# Restore SHA-256-verified official DrivenData GEMS & public USGS/DOE external inputs
# via zero-auth GitHub codeload streams (works both in sandbox and on unrestricted machines).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${ROOT}/.venv/bin/python"
if [ ! -x "${PY}" ]; then
  PY="python3"
fi
exec "${PY}" "${ROOT}/scripts/restore_data.py" "$@"
