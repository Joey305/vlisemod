# Git and deployment hygiene

The repository tracks source code, web assets, the reviewer reproducibility
source package, frozen manifests/inputs/references, and only the small 3EKY
fixture CIF copies. Downloaded PDB corpora, caches, runtime databases,
scientific pipeline outputs, reviewer outputs, and generated release archives
are intentionally local-only.

Before a commit, run:

```bash
python scripts/check_git_hygiene.py --staged
```

The guard rejects forbidden runtime paths and files above 50 MiB unless an
exact reviewed path is supplied with `--allow`. It does not modify the index or
global Git configuration.

Heroku also reads `.slugignore`. It excludes the scientific rebuild and
reviewer-development trees from the web slug; the production Procfile starts
`app:app`, which does not import either tree. The web app can use its configured
RANDY service and remote coordinate retrieval rather than a committed PDB
corpus.
