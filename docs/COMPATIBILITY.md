# BioMate Compatibility Matrix

This document tracks tested versions and known compatibility constraints for BioMate.

## Tested Python Versions

| Python Version | Status | Notes |
|---|---|---|
| 3.12.x | ✅ Tested | Full support |
| 3.13.11 | ✅ Tested (Current) | Full support |
| 3.14.x | ❓ Not tested | Expected to work; pending verification |
| 3.11 and earlier | ❌ Not supported | Uses `itertools.batched` (Python 3.12+) |

**Minimum Requirement:** Python 3.12+ (due to `itertools.batched` usage in sequence generation)

## Core Dependencies

### Updated: NumPy 2.x Compatibility (RESOLVED)

| Package | Tested Version | Status | Notes |
|---|---|---|---|
| numpy | 2.5.0 | ✅ Tested | Migrated from deprecated `numpy.random` to `np.random.Generator` (see ISSUE-002 for details) |

### Tested Dependency Versions

| Package | Tested Version | Status | Purpose |
|---|---|---|---|
| polars | 1.42.0 | ✅ Tested | Data manipulation, CSV/FASTQ parsing |
| altair | 6.2.2 | ✅ Tested | Visualizations (pinned to >=6.2.0,<7.0.0 per ISSUE-003) |
| tornado | 6.5.7 | ✅ Tested | Web interface server |
| regex | 2026.5.9 | ✅ Tested | Pattern matching for index evaluation |
| dnaio | 1.2.4 | ✅ Tested | FASTQ/FASTA file I/O |
| jellyfish | 1.2.1 | ✅ Tested | String distance calculations |

## Known Compatibility Issues

### ISSUE-002: NumPy Random API Deprecation (RESOLVED)

**Status:** ✅ Fixed in commit 8d42be0

- **Problem:** NumPy 2.x deprecated `numpy.random` module (will be removed in 3.x)
- **Solution:** Migrated all random API calls to `np.random.Generator`
- **Migration Pattern:**
  ```python
  # Old (deprecated):
  from numpy import random
  random.choice(items, size=n)
  random.randint(10, 100)
  
  # New (current):
  import numpy as np
  rng = np.random.default_rng()
  rng.choice(items, size=n)
  rng.integers(10, 100)
  ```
- **Affected Modules:** blabber.py (sequence generation)

### ISSUE-003: Altair Version Constraint (RESOLVED)

**Status:** ✅ Fixed in commit f734633

- **Problem:** Loose version spec allowed older untested versions
- **Current Spec:** `altair[all]>=6.2.0,<7.0.0` (pinned to tested version)
- **Affected Modules:** nspector.py (visualization generation)

### ISSUE-004: Docker --break-system-packages Flag (RESOLVED)

**Status:** ✅ Fixed in commit e0f33c0

- **Problem:** Non-standard pip flag masked package conflicts
- **Solution:** Replaced with virtual environment setup
- **Affected Files:** Dockerfile (venv-based installation)

### ISSUE-006: Python 3.12+ Requirement (NO ACTION NEEDED)

**Status:** ✅ Expected behavior

- **Reason:** `itertools.batched` used in blabber.py (line 12) and fastrewind.py
- **Impact:** Only available in Python 3.12+
- **Current Requirement:** pyproject.toml specifies `python = ">=3.12"`
- **Backward Compatibility:** None needed; requirement is intentional

## Breaking Changes & Migration Guide

### NumPy 2.x → 3.x (Future)

When NumPy 3.x is released and `numpy.random` module is removed:

1. ✅ Code is already compatible (migrated to `np.random.Generator`)
2. No changes needed for BioMate when upgrading to NumPy 3.x

### Python 3.12 → 3.13+

No breaking changes detected. Full compatibility confirmed for Python 3.13.11.

## Dependency Update Process

### Safe to Update

The following dependencies can be safely updated:
- **polars**: Major versions tested; minor/patch updates recommended
- **regex**: Stable upstream; updates recommended for new features/fixes
- **dnaio**: Well-maintained; updates recommended
- **jellyfish**: Stable; updates recommended
- **tornado**: Monitor for security updates; apply promptly

### Version-Pinned Dependencies

- **altair**: Pinned to `>=6.2.0,<7.0.0` for stability (requires regression testing on upgrade)
- **numpy**: May require updates when moving between major versions (test `np.random.Generator` API)

### Update Procedure

```bash
# Check for new versions
uv update

# Lock updated versions
uv lock

# Run full test suite
uv run pytest tests/ -v

# If all tests pass, commit new lock file
git add uv.lock
git commit -m "chore: update dependencies"
```

## Testing Recommendations

### Before Upgrading Dependencies

1. Run full test suite: `uv run pytest tests/ -v`
2. Test sequence generation: `uv run biomate blabber --seq-number 1000 --seq-length 100 --output test.fastq`
3. Test visualization: `uv run biomate nspector --input test.fastq --output /tmp`
4. Test sample sheet handling: Verify with real Illumina sample sheets

### CI/CD Environment

BioMate's CI/CD pipeline should test on:
- ✅ Python 3.12 (latest patch)
- ✅ Python 3.13 (latest patch)
- ? Python 3.14 (when available)

## Troubleshooting

### "ImportError: cannot import name 'batched' from 'itertools'"

**Cause:** Running Python < 3.12

**Solution:** Upgrade to Python 3.12 or later:
```bash
python --version  # Check current version
python -m venv --upgrade-deps /path/to/venv  # Update venv
```

### NumPy random API errors

**Old message:** "AttributeError: module 'numpy.random' has no attribute 'randint'"

**Cause:** NumPy 3.x removed deprecated API; code not migrated

**Solution:** Already addressed in BioMate (see ISSUE-002). No action needed.

### Altair visualization not rendering

**Cause:** Version mismatch between altair and dependencies

**Solution:** 
```bash
uv lock --upgrade altair
uv sync
uv run biomate nspector --input test.fastq --output /tmp
```

## Environment Variables

None required. Optional:
- `SSL_CERT` / `SSL_KEY`: For web interface HTTPS (see web_interface.py)

## Reporting Compatibility Issues

If you encounter compatibility problems:

1. Document your environment:
   ```bash
   python --version
   pip list | grep -E "numpy|polars|altair|tornado|regex|dnaio|jellyfish"
   uv lock --frozen  # Show exact versions
   ```

2. Create an issue on GitHub with:
   - Exact error message
   - Steps to reproduce
   - Environment details (Python version, OS, installed versions)
   - Whether issue reproduces on latest main branch

## Version Release Timeline

| Version | Release Date | Status |
|---|---|---|
| 0.3.0 | Current | All known compatibility issues resolved |
| 0.2.x | Earlier | See git history for older tested versions |

---

**Last Updated:** August 2026  
**Maintainer:** BioMate Development Team  
**Next Review:** When adding new dependencies or supporting new Python versions
