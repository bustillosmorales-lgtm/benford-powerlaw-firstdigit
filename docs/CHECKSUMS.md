# Checksums of the raw archives

SHA-256 of the public archives downloaded on 2026-06-27 and used to produce the
paper's numbers. The archives themselves are not redistributed here (see
`QUERIES.md`); verify any re-download against these values.

| SHA-256 | Archive | Source |
|---|---|---|
| `d18c01c615da8f5cd46273b0dacaeee5b1163f4089171da0f84153a5f3c362f6` | `2026q1.zip` | SEC EDGAR Financial Statement Data Sets, 2026Q1 |
| `cc9cd748de4384ab50245bb8d214c989e45649bea6a71ad18a435986ce77a9b2` | `HYG-Database-main.zip` | HYG Database (astronexus) |
| `3d0625177adcac4ec78f12955264bd4a3278e7ba50a5c46ea10332154851e0b3` | `dataverse_files.zip` | MIT Election Lab, Harvard Dataverse VOQCHQ |
| `e8aeeb0cf7001b744371339f2de828228704d10ffb9c52926a012635b6ebc125` | `doi_10_5061_dryad_33f0n__v20200121.zip` | Dryad 10.5061/dryad.33f0n (Laskowski/Pruitt) |
| `f33373eede58e8402ac7e09acd28935db66b388edb021bd0d97f928b6c8492d1` | `osfstorage-archive.zip` | OSF (Gino moral-virtue, Likert) |

Verify with:

```bash
sha256sum <archive>.zip
```
