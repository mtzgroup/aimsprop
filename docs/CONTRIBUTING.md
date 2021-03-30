# Contributing

If you would like to add code to `aimsprop` please follow the steps below.

## Steps to setup package for development

- Pull package from github.

  ```sh
  git clone git@github.com:mtzgroup/aimsprop.git
  cd aimsprop
  ```

- Create a virtual environment and activate it.

  ```sh
  python3 -m venv env
  source ./env/bin/activate
  ```

- Install `flit`. Flit is used to manage package dependencies

  ```sh
  python3 -m pip install flit
  ```

- Source the virtual environment again to ensure you are using your local `flit` install

  ```sh
  source ./env/bin/activate
  ```

- Install the package

  ```sh
  flit install --deps develop --symlink
  ```

- Source one more time to get the newly installed `pytest` cli

  ```sh
  source ./env/bin/activate
  ```

- Make sure tests are passing to test your installation

  ```sh
  pytest
  ```

- Create a new branch for your feature

  ```sh
  git checkout -b feature-{insert-your-great-feature-name-here}
  ```

- Make code contributions and push to GitHub

- Open a Pull Request on GitHub and request a review from one of the package maintainers

## Documentation

- Follow the same steps as above
- Add documentation to the `/docs` folder using [markdown](https://www.markdownguide.org/cheat-sheet) syntax.
- Use `mkdocs` to view your documentation locally. You should be able to view in docs in a browser at [http://127.0.0.1:8000](http://127.0.0.1:8000)

```sh
mkdocs serve
```
