# Workflow Runner

## Overview

This Python project implements a toy parallel workflow runner based on a given JSON workflow specification.

### Approach

The basic approach for this runner is as follows:
- Load DAG workflow from a JSON input (text or file).
- Validate the DAG for issues such as cycles, missing start nodes, type errors, etc. prior to
  running.
- Traverse the DAG using breadth-first search, beginning with a "start" node.
  - Process the node. Here we just print the node's name and mark it visited so that it doesn't get
    processed multiple times.
  - For each node edge, create a task to navigate to the associated child node
    - Use `asyncio.gather` to collect those tasks and run them concurrently.
    - Use `asyncio.sleep` to sleep individual tasks so that other awake tasks can proceed.
    - Once the child node is reached, process the node as described above.
- Run the workflow either directly by importing the python package, or from the Click command line
  interface.

### Workflow Specification

The workflow is specified in the form of a DAG represented in JSON.  The structure consists of
key-value pairs where each key is a node name and the value is an object that describes the node's
properties and its outgoing edges.

#### Key Components
* Node Name: The key for each entry in the JSON object is the name of the node.
* Node Properties: Each node has the following properties:
  * start (optional, boolean): Indicates if the node is the start node of the DAG. There must be
    exactly one start node.
  * edges (object): A dictionary where the keys are the names of the target nodes and the values are
    the delays (in seconds) before visiting the target node.

#### Example
Here is a simple example of a JSON input for the workflow runner:

```json
{
    "A": {"start": true, "edges": {"B": 2, "C": 1}},
    "B": {"edges": {"D": 3}},
    "C": {"edges": {"D": 2}},
    "D": {"edges": {}}
}
```

When run, this workflow should yield:

``` sh
A [elapsed time=0.0]
C [elapsed time=1.0]
B [elapsed time=2.01]
D [elapsed time=3.01]
```

Other examples configurations may be found in `tests/fixtures/` and in `tests/test_main.py`.

### Gross Simplifications

This workflow runner is a toy implementation. It is implemented in pure python and uses the
multi-threading library `asyncio`. While this library does not support true parallelism, it can
switch quickly between threaded asynchronous tasks such that they feel like they are running in
parallel.  This of course assumes that tasks are not CPU intensive but are just slow, like calling
out to an external server for results, or `sleep`). Fortunately that's the situation we're in.

For a real-world workflow runner, considerations such as IO concerns, CPU demand, runtime,
task communication, data size, time requirements, and so on would influence the choice of
parallelization. For CPU-bound parallel processing, we would definitely want to do something truly
parallel using `multiprocessing` or `subprocess` to utilize all available cores. And for any large
compute, we would want to use one of AWS cloud compute services.

Also, in the real-world persistent data stores like PostgreSQL might also be used to track workflows and
tasks. See Airflow and its persistence mechanisms for example.

## Requirements

- Python 3.11+
- Poetry 1.8

## Installation

### Poetry
Use this installation method when developing the package or running tests.

1. Install `poetry` following the instruction [here](https://python-poetry.org/docs/#installation)

2. Install package dependencies:
    ```sh
    poetry install
    ```

### Docker
Or, if you'd prefer to use docker and worry about poetry

1. Install `docker`

2. Build the image

``` sh
docker build -t workflow-runner .
```
Note that the docker image is only set up to run the CLI and doesn't yet accept path inputs.

### Wheel
Finally, to install the project as a python package from the pre-built wheel, run

``` sh
pip install path/to/workflow-runner/dist/workflow_runner-0.1.0-py3-none-any.whl
```

## Usage

You can use the workflow runner either by importing python package into your code or by using
the provided CLI. See source code and associated documentation for python level usage.

### Poetry

#### JSON Input
```sh
poetry run workflow-runner -i '{
    "A": {"start": true, "edges": {"B": 2, "C": 1}},
    "B": {"edges": {"D": 3}},
    "C": {"edges": {"D": 5}},
    "D": {"edges": {}}
}'
```
#### File Input
```sh
poetry run workflow-runner -i tests/fixtures/test_workflow.json
```

#### Help
```sh
poetry run workflow-runner --help
```

### Docker
Note that docker only supports string inputs currently
``` sh
docker run --rm workflow-runner --verbose -i '{
    "A": {"start": true, "edges": {"B": 2, "C": 1}},
    "B": {"edges": {"D": 2}},
    "C": {"edges": {"D": 1}},
    "D": {"edges": {}}
}'
```

## Running Tests
``` sh
poetry run pytest
```
