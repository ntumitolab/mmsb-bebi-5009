# Jupyter book template for Julia Jupyter notebooks

## Jupyter Book

[Jupyter book](https://jupyterbook.org/index.html) builds a website from Markdown and Jupyter Notebook files.

## CI/CD

GitHub actions and GitLab CI/CD are setup to build and publish the website whenever changes are committed.

- Execution results are executed on the fly so you can push notebooks with empty output cells and received the results once the pipeline is completed.
- Execution results are cached so only the edited notebooks would be executed to save CI time.

## Auto update

Periodically updating Julia dependencies and make a PR if notebooks are executed successfully.

- For GitHub: [See the instructions](https://github.com/peter-evans/create-pull-request/blob/main/docs/concepts-guidelines.md#triggering-further-workflow-runs) in `create pull request` action about how to trigger further CI workflows.
- For GitLab: you need a `GIT_PUSH_TOKEN` [CI/CD variable](https://docs.gitlab.com/ee/ci/variables/index.html), which is a PAT with `write_repo` access.
