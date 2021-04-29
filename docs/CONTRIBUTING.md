# Contributing

If you would like to add code to `aimsprop` please follow the steps below.

## Steps to setup package for development

1. Pull package from github.

    ```sh
    git clone git@github.com:mtzgroup/aimsprop.git
    cd aimsprop
    ```

1. Create a virtual environment and activate it.

    ```sh
    python3 -m venv env
    source ./env/bin/activate
    ```

4. Install the latest version of `pip`. Since we will be installing `matplotlib`, having the latest version of `pip` can help avoid issues.

    ```sh
    pip install -U pip
    ```

3. Install `flit`. Flit is used to manage package dependencies

    ```sh
    python3 -m pip install flit
    ```

1. Source the virtual environment again to ensure you are using your local `flit` install

    ```sh
    source ./env/bin/activate
    ```

1. Install the package

    ```sh
    flit install --deps develop --symlink
    ```

1. Source one more time to get the newly installed `pytest` cli

    ```sh
    source ./env/bin/activate
    ```

1. Make sure tests are passing to test your installation

    ```sh
    pytest
    ```

1. Create a new branch for your feature

    ```sh
    git checkout -b feature-{insert-your-great-feature-name-here}
    ```

1. Make code contributions and push to GitHub. See note below for how to add features.

1. Open a Pull Request on GitHub merging your branch into the `develop` branch, and request a review from one of the package maintainers

## How to Add Features to aimsprop

`aimsprop` requires that all new code be well written, reviewed by a package maintainer, documented, and tested. That may sound like a lot, but don't worry--we're here to help you and you'll find this is much easier than it may sound 😀. These steps ensure that `aimsprop` remains usable by other lab members and that the codebase is stable and dependable. Join us in writing great AIMS analysis code!

### Steps to Add a Feature

1. Determine what kind of feature you are hoping to add. Features are usually one of the following:
    1. A property computed on each `Frame` of your `Bundle`. Examples include [compute_bond][aimsprop.geom:compute_bond] or [blur_property][aimsprop.blur:blur_property].
    2. A plotting feature. Note that the current plotting functions [plot_scalar][aimsprop.plot:plot_scalar], [plot_vector][aimsprop.plot:plot_vector], and [plot_population][aimsprop.plot:plot_population] should already be able to plot any new properties you compute.
    3. Creating a subset of the `Bundle` by some value of interest, like a subset by time. [subset_by_t][aimsprop.bundle:Bundle.subset_by_t] is an example.
1. Follow the basic outline below for creating your new feature. Note that all new functions require documentation.

#### Property Computations

New property computation functions take as parameters a `Bundle` object, a `key` that acts as the name for the property, and the relevant parameters to make a calculation. See [compute_bond][aimsprop.geom:compute_bond] code below as an example. The function then performs the desired calculation on each [Frame][aimsprop.bundle:Frame] object and sets the value in the `Frame's` properties dictionary using the `key` passed. See the example in the source code below under [Documenting your Function](#documenting-your-function).

#### New Plotting Function

If you feel the existing plot functions [plot_scalar][aimsprop.plot:plot_scalar], [plot_vector][aimsprop.plot:plot_vector], and [plot_population][aimsprop.plot:plot_population] cannot work for your needs, you can add a new plot generation function to the `plot.py` module.

#### Creating a Subset of the Bundle

To add this feature you would add a new method to the [Bundle][aimsprop.bundle:Bundle] object that takes `self`, and whatever other parameters are required to create your subset, and returns a new `Bundle` object containing the `Frames` of interest. See the `subset_by_x` methods on `Bundle` for examples.

#### Documenting Your Function

Each new function, whether a property computation or a new plotting function, requires proper documentation. The code below demonstrates how to properly document a function.

Documentation should include:

1. An initial "docstring" that, in a single sentence, describes what the function does. Optionally, additional prose can be added underneath this docstring. In the example below this corresponds to the text: "Compute the a bond-angle property for a Bundle (in degrees)."
2. The function signature (the `def my_function_name(...)` part) must contain type annotations that describe the data types passed to your new function. This corresponds to the `bundle: Bundle, key: str, ...` part of the example below.
3. You must describe the parameters and return value of your function using the conventions noted in the example below. This allows for documentation to be generated for your new function and describes to other users what basic values your function requires to operate.

::: aimsprop.geom:compute_angle

#### Formatting Your Code

To assist in creating a clean code format we have created a script for you that will sort your imports, format your code, and notify you of any style violations that you need to fix. To run it, execute the following from the root directory:

```sh
bash ./scripts/format.sh
```

You should run this script frequently as you write your code to help keep your formatting clean.

## Documentation

1. Follow the same steps as above
1. Add documentation to the `/docs` folder using [markdown](https://www.markdownguide.org/cheat-sheet) syntax.
1. Use `mkdocs` to view your documentation locally. You should be able to view in docs in a browser at [http://127.0.0.1:8000](http://127.0.0.1:8000)

```sh
mkdocs serve
```
