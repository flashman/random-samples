import asyncio
import click
import json
from pydantic import BaseModel, Field, model_validator
import time
from typing import Dict, Optional


class Node(BaseModel):
    """
    Represents a node in a Directed Acyclic Graph (DAG).

    Attributes:
        name (str): the name of the node.
        start (bool, optional): Indicates if the node is a start node.
        edges (Dict[str, int]): Dictionary of edges with target nodes and delays.
        visited(bool): Indicates if the node has been processed yet.
    """

    name: str
    start: Optional[bool] = False
    edges: Dict[str, int] = Field(default_factory=dict)
    visited: bool = False


class DAG(BaseModel):
    """
    Represents a Directed Acyclic Graph (DAG) aka the Workflow.

    Attributes:
        nodes (Dict[str, Node]): Dictionary of nodes in the DAG.
    """

    nodes: Dict[str, Node]

    @model_validator(mode="before")
    def check_nodes(cls, values) -> Dict[str, Dict[str, Node]]:
        """
        Check that input nodes represent a single DAG. Specifically, check:
          * valid start node
          * all node keys are present
          * there are no cycles

        Args:
            cls: The class being validated.
            values (Dict[str, Node]): The dictionary of nodes.

        Returns:
            Dict[str, Dict[str, Node]]: The validated dictionary of nodes.

        Raises:
            ValueError: If the inputs are invalid.
        """
        cls._check_start_node(values)
        cls._check_edges(values)
        cls._check_for_cycles(values)

        # Manually set node names.
        for k, data in values.items():
            data["name"] = k

        # Wrap the raw nodes dict.
        return {"nodes": values}

    @staticmethod
    def _check_start_node(values) -> None:
        """
        Check that there is exactly one start node.

        Raises:
            ValueError: If there is not exactly one start node.
        """
        start_nodes = [k for k, v in values.items() if v.get("start", False)]
        if len(start_nodes) != 1:
            raise ValueError("There must be exactly one start node")

    @staticmethod
    def _check_edges(values) -> None:
        """
        Check that all edge targets have corresponding node entries.

        Raises:
            ValueError: If any edge target node is missing.
        """
        node_set = set(values.keys())
        target_set = set(target for v in values.values() for target in v.get("edges", {}).keys())
        if not target_set.issubset(node_set):
            raise ValueError("All edge targets must have corresponding nodes")

    @staticmethod
    def _check_for_cycles(values) -> None:
        """
        Check for cycles in the DAG using depth-first search.

        Raises:
            ValueError: If a cycle is detected in the DAG.
        """
        adjacency_list = {node: data["edges"] for node, data in values.items()}

        visited = set()
        rec_stack = set()

        def dfs(node):
            if node in rec_stack:
                raise ValueError(f"Cycle detected at node: {node}")

            if node in visited:
                return

            visited.add(node)
            rec_stack.add(node)

            for neighbor in adjacency_list[node]:
                dfs(neighbor)
            rec_stack.remove(node)

        for node in adjacency_list:
            if node not in visited:
                dfs(node)

    @property
    def start_node(self) -> "Node":
        """
        Returns the start node.

        Returns:
            Node: The start node.
        """
        for node in self.nodes.values():
            if node.start:
                return node
        raise ValueError("No start node found")


class WorkflowRunner:
    """
    Class to run a workflow asynchronously based on a provided DAG structure.

    Attributes:
        dag (DAG): The Directed Acyclic Graph (DAG) representing the workflow.
        verbose (bool): Setting to control print verbosity.
    """

    def __init__(self, json_input: str, verbose: bool = False) -> None:
        """
        Initializes the WorkflowRunner with the provided JSON input.

        Args:
            json_input (str): JSON string representing the DAG structure.
            verbose (bool): If true, print additional information about the traversal.
        """
        self.dag = DAG.model_validate_json(json_input)
        self.verbose = verbose

    async def visit_target_node_after_delay(self, node: Node, delay: int) -> None:
        """
        Waits for a specified delay and then moves to the target node for processing.

        Args:
            node (Node): The node to visit after delay.
            delay (int): The delay in seconds before visiting the node.
        """
        # Suspend the current task, allowing other tasks to run.
        await asyncio.sleep(delay)
        # Now go ahead and begin processing the target node.
        await self.process_node(node)

    async def process_node(self, node: Node) -> None:
        """
        Processes a node and then schedule concurrent tasks to process children nodes.

        For now, processing just means printing the node name.

        Args:
            node (Node): The node to process.
        """
        # Skip node if we've already been here.
        if node.visited:
            return

        # Perform node processing.
        if self.verbose:
            elapsed_time = round(time.time() - self.start_time, 2)
            print(f"{node.name} [elapsed time={elapsed_time}]")
        else:
            print(node.name)

        # Mark node as visited.
        node.visited = True

        # Schedule concurrent tasks to visit and process target(child) nodes.
        await asyncio.gather(
            *[
                self.visit_target_node_after_delay(self.dag.nodes[target], delay)
                for target, delay in node.edges.items()
            ]
        )

    async def run(self) -> None:
        """
        Starts the asynchronous execution of the workflow starting from the designated start node.
        """
        self.start_time = time.time()
        await self.process_node(self.dag.start_node)


@click.command()
@click.option(
    "-i",
    "--json-input",
    type=str,
    required=True,
    help="JSON string or file path representing the DAG structure.",
)
@click.option("-v", "--verbose", is_flag=True, default=False, help="Enable verbose output.")
def cli(json_input: str, verbose: bool) -> None:
    """
    CLI entry point for running the workflow.

    Args:
        json_input (str): JSON string or file path representing the DAG structure.
        verbose (bool): Enable verbose output.
    """
    try:
        # Check if json_input is a file path or a JSON string.
        if json_input.endswith(".json"):
            with open(json_input, "r") as f:
                json_input = f.read()

        # Initialize and run the workflow.
        runner = WorkflowRunner(json_input, verbose=verbose)
        asyncio.run(runner.run())
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    cli()
