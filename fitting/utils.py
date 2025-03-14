from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Iterator, Sequence

from besmarts.assign.hierarchy_assign_rdkit import (
    smarts_hierarchy_assignment_rdkit as SmartsLabelerRdkit,
)
from besmarts.codecs.codec_rdkit import graph_codec_rdkit as GraphCodecRdkit
from besmarts.core.assignments import graph_assignment as GraphAssignment
from besmarts.core.assignments import graph_db as GraphDb
from besmarts.core.codecs import graph_codec as GraphCodec
from besmarts.core.compute import workqueue_local as LocalWorkqueue
from besmarts.core.graphs import structure as Structure, subgraph_to_structure
from besmarts.core.optimization import optimization_iteration as OptimizationIteration, optimization_step as OptimizationStep
from besmarts.core.perception import perception_model as PerceptionModel
from besmarts.core.trees import tree_node as TreeNode
from besmarts.mechanics.fits import (
    forcefield_optimization_strategy as ForceFieldOptimizationStrategy,
)
from besmarts.mechanics.molecular_models import chemical_system as ChemicalSystem, chemical_system_get_node_hierarchy
from besmarts.mechanics.molecular_models import physical_system as PhysicalSystem
from besmarts.mechanics.smirnoff_models import smirnoff_load
from openff.toolkit import ForceField
from rdkit.Chem import RDKFingerprint


def ff_to_csys(
    force_field: ForceField | str | Path,
    graph_codec: GraphCodec | None = None,
) -> ChemicalSystem:
    """
    Create a chemical system from a SMIRNOFF force field

    Parameters
    ----------
    forcefield
        An OpenFF force field object, or the filename of a SMIRNOFF force field.
    """
    if graph_codec is None:
        graph_codec: GraphCodec = GraphCodecRdkit()

    perception_model = PerceptionModel(
        gcd=graph_codec,
        labeler=SmartsLabelerRdkit(),
    )

    if isinstance(force_field, str) or isinstance(force_field, Path):
        force_field: ForceField = ForceField(force_field)

    with NamedTemporaryFile("wt") as f:
        f.write(force_field.to_string())
        chemical_system = smirnoff_load(
            fname=f.name,
            pcp=perception_model,
        )

    return chemical_system


def data_to_graph_assignment(
    smiles: str,
    data: Sequence[Any],
    graph_codec: GraphCodec | None = None,
) -> GraphAssignment:
    """ """
    if graph_codec is None:
        graph_codec: GraphCodec = GraphCodecRdkit()

    molecule_graph = graph_codec.smiles_decode(smiles)
    selection = {}

    idcs = list(molecule_graph.nodes)
    for i in idcs:
        entry = data[i - 1]

        if (i,) not in selection:
            selection[i,] = []

        selection[i,].append(entry)

    return GraphAssignment(smiles, selection, molecule_graph)

@dataclass
class CandidateGenerator:
    graph_codec: GraphCodec = GraphCodecRdkit()

    @staticmethod
    def get_s0(csys: ChemicalSystem, S: TreeNode) -> Structure:
        hidx = chemical_system_get_node_hierarchy(csys, S)
        assert hidx is not None
        topo = hidx.topology
        return subgraph_to_structure(hidx.subgraphs[S.index], topo)

    def split_ff(self, csys: ChemicalSystem, step: OptimizationStep,):
        from besmarts.core import graphs, configs

        new_candidates = {}
        new_candidates_direct = {}
        direct_success = False
        config: configs.smarts_perception_config = step.pcp

        S = step.cluster
        S0 = self.get_s0(csys, S)

        print(f"Attempting to split {S.name}:")
        s0split = graphs.structure_copy(S0)
        if config.splitter.primitives:
            graphs.graph_set_primitives_atom(s0split, config.splitter.primitives)
            graphs.graph_set_primitives_bond(s0split, config.splitter.primitives)
        print("S0:", self.graph_codec.smarts_encode(S0), "split_space:", self.graph_codec.smarts_encode(s0split))

        if not aa:
            print("No matches.")
            step_tracker[tkey][oper] = strategy.cursor
            return

        print(f"Matched N={len(aa)}")
        seen = set()
        extend_config = config.extender.copy()
        extend_config.depth_max = config.splitter.branch_depth_limit
        extend_config.depth_min = config.splitter.branch_depth_min

        # For each node, I present just.. the chemical objective
        # until I can make a case for IC objective accounting

        for seen_i, i in enumerate(aa, 1):
            g = graphs.graph_to_structure(
                icd.graph_decode(G[i[0]]),
                i[1],
                topo
            )
            graphs.structure_extend(extend_config, [g])
            g = graphs.structure_remove_unselected(g)
            gs = graphs.structure_copy(g)
            if config.splitter.primitives:
                graphs.graph_set_primitives_atom(gs, config.splitter.primitives)
                graphs.graph_set_primitives_bond(gs, config.splitter.primitives)

            if seen_i < 100:
                print(
                    f"{seen_i:06d} {str(i):24s}",
                    # objective.report([x]),
                    gcd.smarts_encode(gs), "<",
                    gcd.smarts_encode(g),
                )
                seen.add(gs)
        print()
        if len(seen) < 2 and len(assn_s) < 100:
            print(f"Skipping {S.name} since all graphs are the same")
            step_tracker[tkey][oper] = strategy.cursor
            return

        if len(seen) < 2 and len(assn_s) < 100:
            print(f"Skipping {S.name} since all graphs are the same")
            step_tracker[tkey][oper] = strategy.cursor
            return

        if graphs.structure_max_depth(S0) > config.splitter.branch_depth_min:
            print("This parameter exceeds current depth. Skipping")
            step_tracker[tkey][oper] = strategy.cursor
            return

        # this is where I generate all candidates
        if step.direct_enable and (
                config.splitter.bit_search_limit < len(assn_s)
            ):
            assn_i = []
            if len(assn_s) < step.direct_limit:
                # or form matches based on unique smarts
                a = []
                extend_config = config.extender.copy()
                extend_config.depth_max = config.splitter.branch_depth_limit
                extend_config.depth_min = config.splitter.branch_depth_min
                for seen_i, (i, x) in enumerate(assn_s.items(), 1):
                    g = graphs.graph_to_structure(
                        icd.graph_decode(G[i[0]]),
                        i[1],
                        topo
                    )
                    graphs.structure_extend(extend_config, [g])
                    g = graphs.structure_remove_unselected(g)
                    a.append(g)

                assn_i = list(range(len(a)))

                pcp = step.pcp.copy()
                pcp.extender = cfg
                print("Direct splitting....")

                ret = splits.split_all_partitions(
                    topo,
                    pcp,
                    a,
                    assn_i,
                    gcd=gcd,
                    maxmoves=0,
                )

                for p_j, (Sj, Sj0, matches, unmatches) in enumerate(ret.value, 0):
                    print(f"Found {p_j+1} {gcd.smarts_encode(Sj)}")

                    edits = 0
                    matches = [
                        y
                        for x, y in enumerate(aa)
                        if x in matches
                    ]
                    unmatches = [
                        y
                        for x, y in enumerate(aa)
                        if x in unmatches
                    ]
                    matched_assn = tuple((assn_i[i] for i in matches))
                    unmatch_assn = tuple((assn_i[i] for i in unmatches))

                    new_candidates_direct[(step.overlap[0], None, p_j, oper)] = (
                        S,
                        graphs.subgraph_as_structure(Sj, topo),
                        step,
                        matches,
                        unmatches,
                        matched_assn,
                        unmatch_assn,
                    )
                if len(new_candidates_direct):
                    direct_success = True
                else:
                    print("Direct found nothing")


        if step.iterative_enable and not direct_success:

            (Q, new_candidates) = union_cache.get((S.index, S.name, strategy.cursor), (None, None))
            # here I need to get the graphs from the gdb
            # where aa refers to the indices
            if Q is None:
                if len(aa) < 100:
                    Q = mapper.union_list_parallel(
                        G, aa, topo,
                        reference=S0,
                        max_depth=graphs.structure_max_depth(S0),
                        icd=icd
                    )
                else:
                    Q = mapper.union_list_distributed(
                        G, aa, topo, wq,
                        reference=S0,
                        max_depth=graphs.structure_max_depth(S0),
                        icd=icd
                    )
                union_cache[(S.index, S.name, strategy.cursor)] = (Q, new_candidates)
            else:
                print(
                    f"{datetime.datetime.now()} Candidates N={len(new_candidates)} retreived from cache for node {S.index}:{S.name}"
                )

            t = datetime.datetime.now()
            print(f"{t} Union is {gcd.smarts_encode(Q)}")

            if new_candidates is None:
                return_matches = config.splitter.return_matches
                config.splitter.return_matches = True

                ret = splits.split_structures_distributed(
                    config.splitter,
                    S0,
                    G,
                    aa,
                    wq,
                    icd,
                    Q=Q,
                )
                config.splitter.return_matches = return_matches

                print(
                    f"{datetime.datetime.now()} Collecting new candidates"
                )
                new_candidates = clusters.clustering_collect_split_candidates_serial(
                    S, ret, step, oper
                )
                for k in new_candidates:
                    v = list(new_candidates[k])
                    v[1] = graphs.subgraph_as_structure(
                        new_candidates[k][1],
                        topo
                    )
                    new_candidates[k] = tuple(v)
                union_cache[(S.index, S.name, strategy.cursor)] = (Q, new_candidates)

        p_j_max = -1
        if candidates:
            p_j_max = max(x[2] for x in candidates) + 1
        for k, v in new_candidates.items():
            k = (k[0], k[1], k[2]+p_j_max, k[3])
            candidates[k] = v
        new_candidates = None

        p_j_max = -1
        if candidates:
            p_j_max = max(x[2] for x in candidates) + 1
        for k, v in new_candidates_direct.items():
            k = (k[0], k[1], k[2]+p_j_max, k[3])
            candidates[k] = v
        new_candidates_direct = None
