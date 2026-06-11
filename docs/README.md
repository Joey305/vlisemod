# V-LiSEMOD Documentation

This `docs/` folder keeps detailed technical, deployment, maintenance, and manuscript-planning material out of the public landing page while preserving a GitHub-friendly map for reviewers and collaborators.

## Documentation Index

| File | Purpose |
|---|---|
| [APP_GUIDE.md](APP_GUIDE.md) | Public-facing explanation of the main app modules and how to interpret each workflow |
| [DATABASE.md](DATABASE.md) | Database layers, key tables, generated assets, and public/private data boundaries |
| [PROTACABILITY.md](PROTACABILITY.md) | Careful interpretation guidance for the PROTACability workflow and outputs |
| [DEPLOYMENT.md](DEPLOYMENT.md) | Local run, environment variables, local-vs-RANDY modes, and deployment caveats |
| [MAINTENANCE.md](MAINTENANCE.md) | Repo hygiene, generated-folder handling, and pre-release checks |
| [DEVELOPER_NOTES.md](DEVELOPER_NOTES.md) | Flask structure, route conventions, feature flags, and developer caveats |
| [MANUSCRIPT_OUTLINE.md](MANUSCRIPT_OUTLINE.md) | Manuscript-ready high-level outline and safe language scaffold |

## Manuscript Planning Folder

The dedicated manuscript-planning material lives under [`docs/manuscript/`](manuscript/):

- [MANUSCRIPT_PLAN.md](manuscript/MANUSCRIPT_PLAN.md)
- [CLAIMS_AND_LIMITATIONS_MATRIX.md](manuscript/CLAIMS_AND_LIMITATIONS_MATRIX.md)
- [FIGURE_AND_TABLE_PLAN.md](manuscript/FIGURE_AND_TABLE_PLAN.md)
- [VALIDATION_AND_REPRODUCIBILITY_PLAN.md](manuscript/VALIDATION_AND_REPRODUCIBILITY_PLAN.md)

## Suggested Reading Order

1. Start with [APP_GUIDE.md](APP_GUIDE.md) for product-level orientation.
2. Read [PROTACABILITY.md](PROTACABILITY.md) before describing degrader-readiness outputs.
3. Use [DATABASE.md](DATABASE.md), [DEPLOYMENT.md](DEPLOYMENT.md), and [MAINTENANCE.md](MAINTENANCE.md) for implementation and operational context.
4. Use [DEVELOPER_NOTES.md](DEVELOPER_NOTES.md) when changing routes, templates, or feature flags.
5. Use [MANUSCRIPT_OUTLINE.md](MANUSCRIPT_OUTLINE.md) and the manuscript folder when drafting publication-facing language.
