
## Contributing

We welcome all contributions!  Please head to [issues](https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects/issues) to either file any bugs/feature requests or find a task you want to assist with.  Make sure to assign yourself the task if you decide to work on it.


### Fork and clone this repository

See the [Github docs](https://help.github.com/articles/fork-a-repo/) for how to make a copy (a fork) of a repository to your own Github account.

Then, [clone the repository](https://help.github.com/articles/cloning-a-repository/) to your local machine so you can begin making changes.

Add this repository as an [upstream remote](https://help.github.com/en/articles/configuring-a-remote-for-a-fork) on your local git repository so that you are able to fetch the latest commits.

On your local machine make sure you have the latest version of the `develop` branch:

```
git checkout develop
git pull upstream develop
```

### Install development dependencies
This will install all the dependencies of the package including the active branch of `aacrgenie`.

```
pip install -r requirements-dev.txt
```


### The development life cycle

`geniesp` follows the standard [git flow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) development strategy.

> To ensure the most fluid development, try not to push to your `develop` or `main` branch.

1. (Assuming you have followed all 4 steps above in the "fork and clone this repository" section). Navigate to your cloned repository on your computer/server.
1. Make sure your `develop` branch is up to date with the `Sage-Bionetworks/GENIE-Sponsored-Projects` `develop` branch.

    ```
    cd {your-github-username}/GENIE-Sponsored-Projects
    git checkout develop
    git pull upstream develop
    ```

1. Create a feature branch which off the `develop` branch. If there is a GitHub/JIRA issue that you are addressing, name the branch after the issue with some more detail (like `{GH|JIRA}-123-add-some-new-feature`).

    ```
    git checkout develop
    git checkout -b GH-123-add-some-new-feature
    ```

1. At this point, you have only created the branch locally, you need to push this to your fork on GitHub.

    ```

    git push --set-upstream origin GH-123-add-some-new-feature
    ```

You should now be able to see the branch on GitHub. Make commits as you deem necessary. It helps to provide useful commit messages - a commit message saying 'Update' is a lot less helpful than saying 'Remove X parameter because it was unused'.

    ```
    git commit changed_file.txt -m "Remove X parameter because it was unused"
    git push
    ```

1. (Make sure you have follow instructions in "Install development dependencies") Once you have made your additions or changes, make sure you write tests and run the test suite.  More information on testing below.

    ```
    pytest ./test
    ```

1. (Make sure you have follow instructions in "Install development dependencies") Make sure to run the auto python code formatter, black.

    ```
    black ./
    ```

1. Once you have completed all the steps above, in Github, create a pull request from the feature branch of your fork to the `develop` branch of Sage-Bionetworks/synapsemonitor.

> *A code maintainer must review and accept your pull request.* A code review ideally happens with both the contributor and the reviewer present, but is not strictly required for contributing. This can be performed remotely (e.g., Zoom, Hangout, or other video or phone conference).

<!-- This package uses [semantic versioning](https://semver.org/) for releasing new versions. The version should be updated on the `develop` branch as changes are reviewed and merged in by a code maintainer. The version for the package is maintained in the [GENIE-Sponsored-Projects/__version__.py](genie/__version__.py) file.  A github release should also occur every time `develop` is pushed into `main` and it should match the version for the package. -->

### Testing

Please add tests for new code. These might include unit tests (to test specific functionality of code that was added to support fixing the bug or feature), integration tests (to test that the feature is usable - e.g., it should have complete the expected behavior as reported in the feature request or bug report), or both.

This package uses [`pytest`](https://pytest.org/en/latest/) to run tests. The test code is located in the [tests](./tests) subdirectory.

Here's how to run the test suite:

```
pytest -vs tests/
```

Tests are also run automatically by Github Actions on any pull request and are required to pass before merging.


### Release Procedure (For Package Maintainers)

Follow gitflow best practices as linked above.

* Always merge all new features into `develop` branch first (unless it is a documentation, readme, or github action patch into `main`)
* After initial features are ready in the `develop` branch, create a `release-X.X` branch to prepare for the release.
    <!-- * update `genie/__version__.py` -->
* Merge `release-X.X` branch into `main`
* Create release tag (`v...`) and include release notes.  Also include any known bugs for each release here.
* Merge `main` back into `develop`

### DockerHub

This repository does not use github actions to push docker images.  By adding the `sagebiodockerhub` github user as an Admin to this GitHub repository, we can configure an automated build in DockerHub.  You can view the builds [here](https://hub.docker.com/r/sagebionetworks/geniesp/builds).  To get admin access to the DockerHub repository, ask Sage IT to be added to the `genieadmin` DockerHub team.
