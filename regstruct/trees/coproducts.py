"""The extraction-contraction coproduct ``δ`` (the negative/renormalization side).

Implements ``D⁻`` of tourist_guide.tex 5636 (= Bruned arXiv:1710.10634 eq. 4.1),
restricted by the ``p₋`` projection to the renormalization coaction

    δ τ = Σ_{φ ⊆ τ,  every component |φ_j|' < 0}  φ  ⊗  τ/^red φ

(plus the unit term ``𝟙₋ ⊗ τ``).  ``φ`` is an extracted sub-forest; ``τ/^red φ``
contracts each component to a red node carrying ``o`` = naive homogeneity of the
extracted component.  For each ``φ`` we sum over

  * a binomial split ``n_φ ≤ n`` of node decorations (coeff ``binom(n, n_φ)``), and
  * boundary-edge decorations ``e'`` on ``∂φ`` (coeff ``1/e'!``), pushed onto the
    extracted nodes (``π e'``) and added to the contracted edges (``e + e'``).

Keeping only divergent components makes the ``e'`` sum finite (raising ``e'``
raises ``|φ|'`` past 0).  Output: a dict ``{(forest, contracted_tree): coeff∈ℚ}``,
``forest`` a canonically-sorted tuple of trees (a product in ``U⁻``).
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import product
from math import comb, factorial

from ..core.homogeneity import Homogeneity, MultiIndex
from .tree import DecoratedTree, red_node, tree

# A coproduct value: {(forest: tuple[DecoratedTree, ...], right: DecoratedTree): coeff}
TensorSum = dict

_E_CAP = 4   # ponytail: safety cap on |e'|_𝔰 per boundary edge; the |φ|<0 filter is the real bound


# --------------------------------------------------------------------------- #
# multi-index helpers
# --------------------------------------------------------------------------- #

def _mi_le(m, n):
    return all(mi <= ni for mi, ni in zip(m, n))


def _mi_binom(n, m):
    r = 1
    for ni, mi in zip(n, m):
        r *= comb(ni, mi)
    return r


def _mi_fact(e):
    r = 1
    for ei in e:
        r *= factorial(ei)
    return r


def _mi_add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def _submultiindices(n):
    """All m with 0 ≤ m ≤ n componentwise."""
    return product(*(range(ni + 1) for ni in n))


def _edge_decs(width, scaling, max_scaled):
    """All multi-indices of length `width` with scaled degree ≤ max_scaled."""
    def rec(i, rem, acc):
        if i == width:
            yield tuple(acc)
            return
        w = scaling.weights[i]
        k = 0
        while w * k <= rem:
            yield from rec(i + 1, rem - w * k, acc + [k])
            k += 1
    yield from rec(0, max_scaled, [])


# --------------------------------------------------------------------------- #
# node-addressed form
# --------------------------------------------------------------------------- #

class _Node:
    __slots__ = ("id", "node_type", "node_dec", "color", "o", "parent", "edge")
    def __init__(self, i, nt, nd, color, o):
        self.id, self.node_type, self.node_dec, self.color, self.o = i, nt, nd, color, o
        self.parent = None      # parent node id
        self.edge = None        # (component, p) of the edge from parent


def _explode(t: DecoratedTree):
    nodes, edges = [], []       # edges: (parent_id, child_id, comp, p)
    def walk(sub, parent_id, edge):
        nid = len(nodes)
        nd = _Node(nid, sub.node_type, sub.node_dec, sub.color, sub.o)
        nodes.append(nd)
        if parent_id is not None:
            nd.parent, nd.edge = parent_id, edge
            edges.append((parent_id, nid, edge[0], edge[1]))
        for (c, p, child) in sub.children:
            walk(child, nid, (c, p))
        return nid
    walk(t, None, None)
    return nodes, edges


def _build(root_id, node_overrides, child_map):
    """Rebuild a DecoratedTree from node-addressed data (recursive)."""
    nt, nd, color, o = node_overrides[root_id]
    children = []
    for (cid, comp, p) in child_map.get(root_id, ()):
        children.append((comp, p, _build(cid, node_overrides, child_map)))
    return tree(nt, nd, children, color=color, o=o)


# --------------------------------------------------------------------------- #
# the coproduct
# --------------------------------------------------------------------------- #

def delta_minus(t: DecoratedTree, sig) -> TensorSum:
    nodes, edges = _explode(t)
    width = sig.width
    n = len(nodes)
    out: TensorSum = defaultdict(Fraction)

    # adjacency
    children_of = defaultdict(list)
    for (a, b, comp, p) in edges:
        children_of[a].append((b, comp, p))

    red_ids = {nd.id for nd in nodes if nd.color == "red"}

    for mask in range(1 << n):
        phi = {i for i in range(n) if mask & (1 << i)}
        if not red_ids <= phi:        # φ must contain all red nodes
            continue
        if not phi:                   # unit term 𝟙₋ ⊗ τ
            out[((), t)] += Fraction(1)
            continue

        # connected components of the induced subgraph on φ
        comp_id = {}
        for i in sorted(phi):
            par = nodes[i].parent
            comp_id[i] = comp_id[par] if (par in phi) else i  # parent processed first (pre-order ids)
        comps = defaultdict(list)
        for i in phi:
            comps[comp_id[i]].append(i)

        # boundary edges ∂φ: (x∈φ) → (y∉φ)
        boundary = [(a, b, comp, p) for (a, b, comp, p) in edges if a in phi and b not in phi]

        # node-decoration split n_φ ≤ n on φ-nodes
        dec_choices = {i: list(_submultiindices(nodes[i].node_dec)) for i in phi}
        # boundary-edge decoration e' choices
        edge_choices = list(_edge_decs(width, sig.scaling, _E_CAP))

        for nphi in product(*(dec_choices[i] for i in sorted(phi))):
            nphi_map = dict(zip(sorted(phi), nphi))
            binom = 1
            for i in phi:
                binom *= _mi_binom(nodes[i].node_dec, nphi_map[i])
            if binom == 0:
                continue
            for eprime in product(edge_choices, repeat=len(boundary)):
                ep_map = {k: eprime[k] for k in range(len(boundary))}
                efact = 1
                for ev in eprime:
                    efact *= _mi_fact(ev)
                coeff = Fraction(binom, efact)

                term = _assemble(t, nodes, children_of, phi, comps, comp_id,
                                 boundary, nphi_map, ep_map, sig, width)
                if term is None:
                    continue
                forest, right = term
                # p₋: keep only divergent extracted components
                if any(not c.extended_homogeneity(sig).is_negative() for c in forest):
                    continue
                out[(forest, right)] += coeff

    return {k: v for k, v in out.items() if v != 0}


def _assemble(t, nodes, children_of, phi, comps, comp_id, boundary,
              nphi_map, ep_map, sig, width):
    zero = (0,) * width
    # π e': boundary-edge decoration pushed onto its φ-endpoint
    pi_e = defaultdict(lambda: zero)
    edge_extra = {}             # boundary edge index -> e' (added to the contracted edge)
    for idx, (a, b, comp, p) in enumerate(boundary):
        ev = ep_map[idx]
        pi_e[a] = _mi_add(pi_e[a], ev)
        edge_extra[(a, b)] = ev

    # ---- left forest: each φ-component as a tree (rooted at component top) ----
    forest = []
    for top, members in comps.items():
        member_set = set(members)
        ov = {}
        cmap = defaultdict(list)
        for i in members:
            nd = nodes[i]
            ndec = _mi_add(nphi_map[i], pi_e[i])
            ov[i] = (nd.node_type, ndec, nd.color, nd.o)
            for (cid, comp, p) in children_of[i]:
                if cid in member_set:           # edge interior to the component
                    cmap[i].append((cid, comp, p))
        forest.append(_build(top, ov, cmap))
    forest = tuple(sorted(forest, key=lambda x: x._sortkey()))

    # ---- right tree: contract each component to a red node, keep non-φ nodes ----
    ov = {}
    for nd in nodes:
        if nd.id not in phi:
            ov[nd.id] = (nd.node_type, nd.node_dec, nd.color, nd.o)
    for top, members in comps.items():
        # contracted node decoration  [n − n_φ]_φ  (summed over the component)
        ndec = zero
        o_acc = Homogeneity(0)
        for i in members:
            nd = nodes[i]
            rem = tuple(a - b for a, b in zip(nd.node_dec, nphi_map[i]))
            ndec = _mi_add(ndec, rem)
            if nd.color == "red":
                o_acc = o_acc + nd.o
        # o gains the naive homogeneity of the extracted component
        o_acc = o_acc + _component_naive_homog(members, nodes, nphi_map, pi_e, children_of, sig)
        ov[top] = ("bullet", ndec, "red", o_acc)

    # rebuild the right tree: representative id of a φ-node is its component top
    rep = {i: (comp_id[i] if i in phi else i) for i in range(len(nodes))}
    cmap = defaultdict(list)
    for (a, b, comp, p) in [(e[0], e[1], e[2], e[3]) for e in _edges_of(children_of)]:
        ra, rb = rep[a], rep[b]
        if ra == rb:                       # interior to a component → not in right
            continue
        extra = edge_extra.get((a, b), zero)
        cmap[ra].append((rb, comp, _mi_add(p, extra)))
    root_rep = rep[0]
    right = _build(root_rep, ov, cmap)
    return forest, right


def _edges_of(children_of):
    for a, lst in children_of.items():
        for (b, comp, p) in lst:
            yield (a, b, comp, p)


def _component_naive_homog(members, nodes, nphi_map, pi_e, children_of, sig):
    """Naive homogeneity |component|' with decorations n_φ + π e' and interior edges."""
    member_set = set(members)
    h = Homogeneity(0)
    for i in members:
        nd = nodes[i]
        base = Homogeneity(0) if nd.color == "red" else sig.node_homogeneity(nd.node_type)
        ndec = _mi_add(nphi_map[i], pi_e[i])
        h = h + base + Homogeneity(sig.scaled(ndec))
        for (cid, comp, p) in children_of[i]:
            if cid in member_set:
                h = h + sig.edge_gain(comp, p)
    return h


# --------------------------------------------------------------------------- #
# the group coproduct δ⁻ : U⁻ → U⁻ ⊗ U⁻  (right factor projected by p₋)
# --------------------------------------------------------------------------- #

def _p_minus(t: DecoratedTree, sig):
    """p₋ on a single tree: keep if divergent; ● and red ●^{0,α} → 𝟙₋; else 0."""
    if t.homogeneity(sig).is_negative():
        return ("keep", t)
    if t.nodes() == 1 and not any(t.node_dec) and (t.color == "red" or t.node_type == "bullet"):
        return ("unit",)
    return ("zero",)


def _sort_forest(f):
    return tuple(sorted(f, key=lambda x: x._sortkey()))


def delta_minus_group(t: DecoratedTree, sig):
    """δ⁻ on a single generator tree → {(left_forest, right_forest): coeff}.

    ponytail: WIP — this naive per-step p₋ collapse does not implement the BHZ
    Hopf-algebra *quotient* consistently; coassociativity fails on the two-edge
    gKPZ tree (the extended-decoration re-extraction subtlety, BHZ Remark 5.38).
    Single-application δ (`delta_minus`) is validated; this needs the quotient
    done properly (work in the extended space T̂⁻, then quotient).
    """
    out = defaultdict(Fraction)
    for (left, right), coeff in delta_minus(t, sig).items():
        pm = _p_minus(right, sig)
        if pm[0] == "zero":
            continue
        rf = (right,) if pm[0] == "keep" else ()
        out[(left, rf)] += coeff
    return {k: v for k, v in out.items() if v != 0}


def _delta_group_forest(forest, sig):
    """δ⁻ extended to a forest (algebra morphism): δ⁻(∏ τ_i) = ∏ δ⁻(τ_i)."""
    acc = {((), ()): Fraction(1)}
    for comp in forest:
        d = delta_minus_group(comp, sig)
        new = defaultdict(Fraction)
        for (l1, r1), c1 in acc.items():
            for (l2, r2), c2 in d.items():
                new[(_sort_forest(l1 + l2), _sort_forest(r1 + r2))] += c1 * c2
        acc = dict(new)
    return acc


def coassoc_lhs(t, sig):
    """(δ⁻ ⊗ id) δ⁻ τ  →  {(forest, forest, forest): coeff}."""
    out = defaultdict(Fraction)
    for (A, B), c in delta_minus_group(t, sig).items():
        for (A1, A2), c2 in _delta_group_forest(A, sig).items():
            out[(A1, A2, B)] += c * c2
    return {k: v for k, v in out.items() if v != 0}


def coassoc_rhs(t, sig):
    """(id ⊗ δ⁻) δ⁻ τ  →  {(forest, forest, forest): coeff}."""
    out = defaultdict(Fraction)
    for (A, B), c in delta_minus_group(t, sig).items():
        for (B1, B2), c2 in _delta_group_forest(B, sig).items():
            out[(A, B1, B2)] += c * c2
    return {k: v for k, v in out.items() if v != 0}


# --------------------------------------------------------------------------- #
# the negative twisted antipode S'₋ : U⁻ → ℝ[U]   (tourist_guide.tex 5034)
# --------------------------------------------------------------------------- #

def twisted_antipode(t: DecoratedTree, sig, memo=None):
    """S'₋(τ) as a forest-sum {forest: coeff} in ℝ[U]; algebra morphism, recursive.

        S'₋(τ) = −M₋(S'₋⊗Id)(δτ − τ⊗●^{0,|τ|}).
    """
    if memo is None:
        memo = {}
    key = t._sortkey()
    if key in memo:
        return memo[key]
    d = delta_minus(t, sig)
    top = ((t,), red_node(t.extended_homogeneity(sig), width=sig.width))
    result = defaultdict(Fraction)
    for (left, right), coeff in d.items():
        if (left, right) == top:                 # subtract τ ⊗ ●^{0,|τ|}
            continue
        for f, c2 in _antipode_forest(left, sig, memo).items():
            result[_sort_forest(f + (right,))] += coeff * c2
    result = {f: -v for f, v in result.items() if v != 0}
    memo[key] = result
    return result


def _antipode_forest(forest, sig, memo):
    """S'₋ extended to a forest (algebra morphism)."""
    acc = {(): Fraction(1)}
    for comp in forest:
        sc = twisted_antipode(comp, sig, memo)
        new = defaultdict(Fraction)
        for f1, c1 in acc.items():
            for f2, c2 in sc.items():
                new[_sort_forest(f1 + f2)] += c1 * c2
        acc = dict(new)
    return acc
