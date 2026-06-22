# Code Quality Audit — `counterterms/`

**Date:** 2026-06-22
**Scope:** All of `counterterms/` **except** `render/` and `tests/` (per request). Code quality only —
**not** mathematical correctness, **not** test coverage. Read-only: no code was changed.
**Method:** four parallel readers (core+structures, equation, trees, renorm) plus cross-cutting
static checks (`pyflakes`, `radon`, grep). Each finding below was spot-verified against the source.

## Verdict

Healthy research codebase. Strong, paper-cited docstrings; a clean layered design
(`homogeneity → signature/symbol → hopf → structures → api`); honest `[SOCKET]` /
`NotImplementedError` markers at deferred boundaries; every module has `from __future__ import
annotations`; no `TODO`/`FIXME`/`XXX` left lying around. `homogeneity.py`, `hopf.py`,
`signature.py`, `rule.py`, and `nonlinearity.py` are the cleanest.

Two systemic issues and a handful of localized ones are worth addressing. Nothing is a
correctness defect (not in scope, but none surfaced incidentally either).

---

## Systemic issue 1 — incomplete type annotations (recurring, project-mandate violation)

CLAUDE.md requires **full, precise** type hints on every function/method param and return, and on
dataclass fields, preferring a precise domain type over bare containers. The single highest-volume
gap across the package is unparameterized `dict` / `list` / `tuple` / `set` / `object` where the real
shape is known. The values are almost always one of: a tensor-sum `dict[tuple[...], Fraction]`, a
free-module `dict[DecoratedTree, ...]`, `list[Counterterm]`, `list[Unknown]`, `list[DecoratedTree]`,
or `sympy.Expr`.

Worst offenders:

| File | Lines | What's bare |
|------|-------|-------------|
| `trees/coproducts.py` | `36` | `TensorSum = dict` — a type alias that throws away all shape info; applied to most public fns yet sibling fns use bare `-> dict` for the same thing (inconsistent). |
| `structures.py` | `39,40,108,177` fields; `42,45,48,110,117,123,127,131,148,180,188,202,206,211` returns | pervasive bare `dict`/`tuple`/`list`/`set`/`object`. |
| `renorm/equation.py` | `47–50,57,61,70` | `base: dict`, `unknowns: list`, `per_component: dict`, `all_trees: list`, `counterterms()->list`, `_display_subs()->dict`. `Unknown`/`Counterterm`/`DecoratedTree` not in the `TYPE_CHECKING` block. |
| `equation/dsl.py` | `46–78` (instance attrs), `83–84`, `129`, `140`, `179`, `181` | `Unknown`/`Noise`/`Parabolic` **instance attributes unannotated**; `SPDE.equations/noises: list`; `build_context -> tuple[Signature, dict, list]`; `base` values typed `object` though always `sympy.Expr`. |
| `core/hopf.py` | `30,48–51,64,84` | generic `Callable[[Hashable], dict]` / bare `dict` returns. |
| `renorm/scheme.py` | `43,58,71,121` | `Iterable`, `tuple`, `list`, `dict` unparameterized. |
| `renorm/nonlinearity.py` | `47` | `base: dict`. |
| `trees/tree.py` | `42,76` | `_sortkey()->tuple`, `canonical_key()->tuple` (the load-bearing key — an alias would help). |
| `equation/generate.py` | `53,108` | `pool: dict`, `counts: Counter`. |

This is mechanical to fix and would meaningfully improve readability of the tensor-sum/free-module
data flow, which is currently all spelled out in docstrings rather than types.

## Systemic issue 2 — complexity / duplication hotspots

Three functions carry their complexity heavily (radon cyclomatic complexity in parens):

- **`trees/coproducts.py:154 delta_minus` (CC 50)** and **`:462 delta_plus` (CC 41)** — each ~100
  lines, 4–5 levels of nesting, and they share substantial **copy-pasted** machinery:
  - the `binom` computation: `211–215` ≡ `506–508` (character-for-character).
  - `efact` from `eprime`: `234–236` ≡ `523–525` (identical).
  - `pi_e` / `edge_extra` construction from the boundary: `266–272` ≈ `527–531`.
  - both `_explode`, build `children_of`, loop `for mask in range(1<<n)`, compute a boundary.

  Extracting 3–4 small shared helpers (`_binom_over`, `_efact`, `_push_eprime`) and pulling the
  self-contained `p₋` projection (`243–255`) and `e′`-prune (`219–232`) out of `delta_minus` into
  named helpers would cut the nesting without touching the math. **High** for readability;
  `delta_minus` is the hardest function in the package to follow.

- **`equation/dsl.py:129 build_context` (CC 31)** — does parse, per-equation noise/`g` split, four
  distinct scope checks, and the structural-rule `allowed` construction all inline. Splitting the
  scope checks and the `allowed`-builder into helpers would help. **Med.**

- **`structures.py:56 bhz_character` vs `:66 canonical_character`** — near-identical copy-paste
  (verified): same loop over `twisted_antipode(t)`, same `Rational`, same forest product, same
  `expand`. The only difference is `canonical_character` zeroes odd-noise forests. `bhz_character`
  *is* `canonical_character` with the parity check disabled — collapse to one with a flag. **High**
  (clearest single refactor target).

---

## Per-file findings

### `core/homogeneity.py` — clean
- **Low** `36,39,48` redundant string forward-refs (`"Homogeneity"`) despite `__future__` annotations.
- **Low** defines `__lt__`/`__eq__` but callers sort via `key=lambda h: h._key()` (structures.py
  `120,125`) rather than trusting `<` — partial comparison surface; `functools.total_ordering` would
  make intent explicit. `_key()` is "private" yet part of the public sort contract across modules.
- **Low** `60` `__str__` packs three nested conditionals into one f-string — densest line in the file.

### `core/jets.py` — fragile string encoding (latent robustness risk)
- **Med** `19,jet_parts` decode SymPy symbol names (`u{comp}_{i}_{j}`) by string parsing with no
  validation — a component index ≥ 10 or a malformed name silently mis-parses. It's a trust boundary
  with no guard, and the name grammar is documented nowhere. A one-line `# u{comp}_{k0}_{k1}_…`
  comment + a guard would remove the ambiguity. `is_jet` checks `s.name[1].isdigit()` while
  `jet_parts` assumes the whole first segment is the int — consistent only for single-digit comps.

### `core/signature.py` — clean
- **Med** `40–43 node_homogeneity` comment "● (and red) nodes have homogeneity 0" is misleading:
  callers special-case `red` *before* calling this; the method handles red only by accident. Stale.
- **Low** `Signature` is a mutable `@dataclass` (unlike `Homogeneity`/`Scaling`) yet is the shared
  context passed everywhere and treated as immutable — `frozen=True` would enforce that.
- **Low** `23` unnecessary quoted forward-ref; `49–50 scaled` is a thin duplicate wrapper over
  `self.scaling.scaled`.

### `core/symbol.py` — speculative single-impl Protocol
- **Med (over-engineering)** single-implementation `Symbol` Protocol; `DecoratedTree` is the only
  impl, the multi-index basis is explicitly "not built". Lightweight (two methods) and honestly
  marked `[SOCKET]`, so cost is low — flagging per brief. Nothing checks `DecoratedTree` actually
  satisfies it (`@runtime_checkable` but no `isinstance` / test), so the Protocol is decorative today.
- **Low** `14` `typing.Hashable` is the deprecated alias of `collections.abc.Hashable` (same in
  `hopf.py:24`).

### `core/hopf.py` — strong, generic
- **Med** loose generic typing on the structural maps (see systemic issue 1).
- **Low** `38` `sum((...), 0)` seeds with int `0` over symbolic values (works via `__radd__`, but
  mixes int into a symbolic accumulation; same at `structures.py:209`). `comodule_action` seeds
  `defaultdict(int)` (`90`) while `antipode` seeds `defaultdict(Fraction)` (`69`) — two zero
  conventions in one file. `74` `Fraction(cl)` re-wrap is redundant (recursion already returns
  `Fraction`).

### `structures.py` — see systemic issues; plus:
- **Med** `202–204 convolve` method shadows the imported `convolve` from `core.hopf` (`23`); and
  `antipode` is imported top-level (`23`) **and** re-imported locally inside `structure_antipode`
  (`138`) — redundant duplicate import + readability trap.
- **Low** non-parallel naming across the two classes: `coproduct` (is actually `δ⁻`), `coaction`,
  `recentering`, `structure_coproduct` — requires reading docstrings to disambiguate.
- **Low** `27` imports the private `_delta_group_forest` from `coproducts` (cross-module private).
- **Low** `209 inverse` rebuilds the antipode (with its memo) on every call — no reuse. Acceptable
  for a precision-over-speed tool.
- **Positive** `admissible()` raising a clear, cited `NotImplementedError` at a deferred boundary.

### `__init__.py` — clean public surface
- **Low** `4` module docstring still says "Phase 1: …" but `__all__` exports Phase 3/4 entry points
  (`build_regularity_structure`, `build_renormalization_group`, `daprato_lift`) — **stale**.
- **Low** `renormalize` wrapper has no docstring while the other three wrappers do; `jet` exported at
  top level is an odd granularity (everything else is a DSL noun or pipeline verb).

### `api.py` — short and clear
- **Low** `21,25` accumulators `per_component`/`contribs` untyped (locals, so minor).
- **Low** `29` `if ed != 0` relies on SymPy truthiness; `ed.is_zero` / compare to `sympy.Integer(0)`
  matches the SymPy-aware style used elsewhere and is more explicit.

### `equation/dsl.py` — see systemic issues; plus:
- **Med (smell)** `81 @dataclass` on `SPDE` with a **fully hand-written `__init__`** (`86–97`) that
  overrides the generated one — the decorator buys nothing but a `__repr__`/`__eq__` over two bare
  fields, and misleads the reader. Either drop `@dataclass` or use the generated init. *(verified)*
- **Med** `134` local `width` assigned but never used *(confirmed by pyflakes)*.
- **Low** `74` `import warnings` inside `Parabolic.__init__` instead of module top.
- **Low** `146` re-`expand`s `rhs` inside the noise loop though it was expanded at `143`.
- **Low** `10` docstring scope `β₀ ∈ (−2,0)` is stale/narrow — subcriticality is rule-based and the
  da Prato lift extends below −2.

### `equation/generate.py` — excellent docstrings
- **Med** `70` `assert guard < 50` is the only non-termination backstop — asserts are stripped under
  `python -O`, so it vanishes in optimized runs. For a research backstop, prefer an explicit `raise`.
- **Low** `53 pool: dict` used as an ordered set with key==value — intentional but worth a one-line
  note. `104 _emit(bound=_BOUND)` default is dead (always called with `bound` passed). `_MAX_NODE_SCALED
  = 2` cap (`32`) is undocumented vs. the `|n|_𝔰 ≤ 1` scope.

### `equation/rule.py` — cleanest file
- **Med** `19` `from fractions import Fraction` imported but unused *(confirmed by pyflakes)*.
- **Low** `28–29` early `return` treats a noiseless signature as subcritical (correct, but the intent
  is implicit — one comment would help). Excellent raising contract + error message otherwise.

### `equation/daprato.py` — good scope-rejection error handling
- **Med** `54 X = sympy.Symbol("_X_")` and `66 f"X{k}"` Wick-noise names can silently collide with a
  user-named field/noise; no guard. Low probability in scope, but un-namespaced magic symbols.
- **Low** `64–65` rational round-trip via `.numerator/.denominator` is unnecessary — `sympy.Rational`
  accepts a `Fraction` directly. `59 deg` guard is mildly redundant.

### `trees/tree.py` — clean, well-cited
- **Low** `40 o: Homogeneity = field(default_factory=lambda: _ZERO)` uses a factory for an immutable
  singleton while `82 tree(... o=_ZERO)` uses the bare default directly — inconsistent treatment of
  the same default in one file. `29 Edge` string forward-ref is actually load-bearing in an alias
  assignment (worth a one-word comment so it doesn't read as a mistake).

### `trees/coproducts.py` — see systemic issues; plus:
- **Med** `90–95 _Node` is a hand-rolled `__slots__` class with a tuple-unpack `__init__`;
  `@dataclass(slots=True)` would be shorter, typed per-field, and consistent with `DecoratedTree`.
- **Med** `312` `[(e[0],e[1],e[2],e[3]) for e in _edges_of(...)]` is a no-op repack of 4-tuples
  `_edges_of` already yields — dead/pointless. `233 ep_map` rebuilds `dict(enumerate(eprime))` for
  no gain (a list/tuple would do). `480–490 _subext` is a memoized recursion that's then eagerly
  forced for every node — plainer than a post-order fill needs.
- **Low** `38 _E_CAP = 4` and its `# ponytail:` comment ("the |φ|<0 filter is the real bound") now
  mildly contradict the later e′-prune blocks (`219–232`, `510–521`) which say *they* are the real
  bound "here made effective" — `_E_CAP` is effectively a redundant safety cap. `264 _assemble`
  parameter is named `boundary` but `delta_minus` passes `recenter` (boundary + between edges) — the
  name is now stale at the callee. `defaultdict(lambda: (0,)*width)` vs `defaultdict(lambda: zero)`
  spelled two ways (`267,527,547`). `555 BLUE = -1` in-band sentinel id.
- **Note (not verified, math out of scope)** `225–231` the e′-prune computes each edge's allowance
  independently against the un-recentered base — a deliberately *loose* (superset) prune, safe
  because the exact `p₋` filter at `243–255` still runs. The independence assumption is subtle and
  undocumented as such — worth a one-line note.

### `renorm/nonlinearity.py` — cleanest of renorm
- **Low** `32` `for s in list(expr.free_symbols)` — the `list(...)` wrapper is unnecessary (set isn't
  mutated). `55 elem_diff` calls `sympy.expand` at every recursion level — redundant on large trees
  (perf smell, not a bug). `_D` has a docstring, `_Dn` doesn't — asymmetric.

### `renorm/equation.py` — see systemic issue 1; plus:
- **Med (duplication)** two parallel jet→display mechanisms: `_display_subs`/`_all_display_subs`
  (`61–74`, → `sympy.Derivative`) and `_pretty` (`86–100`, → compact subscripts). Both walk jets and
  rebuild substitution dicts for different output targets — a latent drift risk.
- **Low** `31 Counterterm.elem_diff` field shadows the module name `elem_diff` from `nonlinearity.py`
  (no real collision; confusing). `58` docstring "Counterterms of component 0 (the scalar case)" is
  accurate but the property name `counterterms` reads as "all". Lazy imports inside methods (`88,112,
  119,130`) are a deliberate cycle-avoidance pattern — fine, worth one comment saying so.

### `renorm/scheme.py` — over-engineering hotspot (but intentional phasing)
- **Med→High (judgment call)** the file ships a full scheme/character/evaluator **injection seam**
  ahead of any real consumer: `RenormalizationScheme` Protocol (`130`), `IntegralEvaluator` Protocol
  (`172`), `UnbuiltEvaluator` (`181`, only ever raises), and `BPHZ.numeric_character` (`164`, can
  only raise). This is exactly the kind of speculative infrastructure ponytail/YAGNI would defer.
  **However**, CLAUDE.md/ROADMAP explicitly list "canonical BPHZ *values* (Wick/integrals)" as the
  next planned track (Phase 4 / Track B2) and the code is *honest* about it via `[SOCKET]` comments.
  So this is intentional, documented phasing rather than accidental bloat — but it is the clearest
  place to delete/defer code if the team decides to apply YAGNI strictly until B2 actually starts.
- **Med** `23` imports the private `_explode` from `..trees.coproducts` across a package boundary —
  should be public if it's a real cross-module dependency.
- **Low** `63 has_odd_noise` and the inline `len(noise_ids) % 2 == 1` in `expectation` (`89`)
  duplicate the parity check + an `_explode` pass. `19 Iterable` from `typing` (deprecated location).

---

## Suggested priority order (if/when acting on this — no changes made here)

1. **Collapse `bhz_character`/`canonical_character`** (structures.py) — clearest single win.
2. **Remove the two dead/unused items** — `rule.py:19` unused `Fraction` import, `dsl.py:134` unused
   `width`; and the no-op repack at `coproducts.py:312`.
3. **Fix `SPDE`'s `@dataclass` + overridden `__init__`** (dsl.py:81) — drop the decorator or the init.
4. **Type-annotation pass** (systemic issue 1) — start with `TensorSum = dict` and the dataclass
   fields/returns in `structures.py`, `renorm/equation.py`, `dsl.py`.
5. **Decompose `delta_minus`/`delta_plus`** into shared helpers (coproducts.py) — biggest readability
   payoff, most care needed (math-adjacent).
6. **Decide on `scheme.py`'s unbuilt seam** — keep as documented phasing, or trim to `FreeConstants`
   + `BPHZ.character` until Track B2 lands.
7. Doc fixes: stale `__init__.py` Phase-1 docstring, `signature.node_homogeneity` "(and red)"
   comment, `dsl.py` `β₀ ∈ (−2,0)` scope line, `_assemble`'s stale `boundary` param name.
