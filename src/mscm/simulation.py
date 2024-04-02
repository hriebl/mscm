"""
TODO.
"""

import jax.numpy as jnp
import jax.random
import liesel.model as lsl
import networkx as nx


def _build_simulation_graph(nodes):
    edges = []

    for node in nodes:
        for _input in node.all_input_nodes():
            if isinstance(node, lsl.Dist) and _input is node.at:
                edges.append((node, _input))
            else:
                edges.append((_input, node))

    graph = nx.DiGraph(edges)
    graph.add_nodes_from(nodes)
    return graph


def _simulate_model(model, seed, skip=()):
    graph = _build_simulation_graph(model.nodes.values())
    nodes = list(nx.topological_sort(graph))

    dists = [node for node in nodes if isinstance(node, lsl.Dist)]
    seeds = jax.random.split(seed, len(dists))

    for dist, seed in zip(dists, seeds):
        if dist.name in skip:
            continue

        tfp_dist = dist.init_dist()
        event_shape = tfp_dist.event_shape
        batch_shape = tfp_dist.batch_shape
        value_shape = jnp.asarray(dist.at.value).shape
        sample_index = len(value_shape) - len(batch_shape) - len(event_shape)
        sample_shape = value_shape[:sample_index]

        if isinstance(dist.at, lsl.TransientIdentity):
            dist.at.inputs[0].value = tfp_dist.sample(sample_shape, seed)
        else:
            dist.at.value = tfp_dist.sample(sample_shape, seed)


def simulate_model(model, seed, skip=()):
    seeds = jax.random.split(jax.random.PRNGKey(seed))

    _simulate_model(model, seeds[0], skip)
    index = jnp.sum(model.vars["z"].value, axis=1) < 0.1
    model.vars["n"].value = model.vars["n"].value.at[index].set(0.0)
    skip = [name for name in model.nodes if name != "y_log_prob"]
    _simulate_model(model, seeds[1], skip)
