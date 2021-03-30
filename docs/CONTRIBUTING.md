# Contributing

If you would like to add code to `aimsprop` please follow the steps below.

## Steps to setup package for development

1. Pull package from github.

    ```sh
    git clone git@github.com:mtzgroup/aimsprop.git
    cd aimsprop
    ```

2. Create a virtual environment and activate it.

    ```sh
    python3 -m venv env
    source ./env/bin/activate
    ```

3. Install `flit`. Flit is used to manage package dependencies

    ```sh
    python3 -m pip install flit
    ```

4. Source the virtual environment again to ensure you are using your local `flit` install

    ```sh
    source ./env/bin/activate
    ```

5. Install the package

    ```sh
    flit install --deps develop --symlink
    ```

6. Source one more time to get the newly installed `pytest` cli

    ```sh
    source ./env/bin/activate
    ```

7. Make sure tests are passing to test your installation

    ```sh
    pytest
    ```

8. Create a new branch for your feature

    ```sh
    git checkout -b feature-{insert-your-great-feature-name-here}
    ```

9. Make code contributions and push to GitHub
10. Open a Pull Request on GitHub and request a review from one of the package maintainers

## Documentation

1. Follow the same steps as above
2. Add documentation to the `/docs` folder using [markdown](https://www.markdownguide.org/cheat-sheet) syntax.
3. Use `mkdocs` to view your documentation locally. You should be able to view in docs in a browser at [http://127.0.0.1:8000](http://127.0.0.1:8000)

```sh
mkdocs serve
```
