Title: Combinatorics of the cosmohedron

URL Source: https://arxiv.org/html/2603.03425v3

Published Time: Mon, 24 Aug 2026 21:34:11 GMT

Markdown Content:
Federico Ardila-Mantilla Nima Arkani-Hamed Note:Institute for Advanced Study, Princeton; arkani@ias.edu.Carolina Figueiredo Note:Jadwin Hall, Princeton University; cfigueiredo@princeton.edu. Francisco Vazão Note:Institute for Advanced Study, Princeton; fvvazao@ias.edu.

###### Abstract

The cosmohedron was recently proposed as a polytope underlying the cosmological wavefunction for \text{Tr}(\Phi^{3}) theory. Its faces were conjectured to be in bijection with Matryoshkas, which are obtained from a subdivision of a polygon by sequentially wrapping groups of polygons into larger polygons. In this paper, we prove the correctness of this construction and elucidate its combinatorial structure. Cosmohedra generalize to a wider class of \mathcal{X} in Y polytopes, where we chisel a polytope from the family \mathcal{X} at each vertex of a polytope Y. We sketch a new application of these chiseled polytopes to the physics of ultraviolet divergences in loop-integrated Feynman amplitudes.

## 1 Introduction

It has been a long-standing question in cosmology to find a geometric object that encodes the cosmological wavefunction of Tr(\phi^{3}) theory. The combinatorics of this wavefunction is governed by certain nesting subdivisions of a polygon that we call _Matryoshkas_.

Arkani-Hamed, Figueiredo, and Vazão [[AHFV25](https://arxiv.org/html/2603.03425v3#bib.bibx8)] proposed a construction of the _cosmohedron_ as a solution to this cosmological question. The goal of this article is to prove the correctness of that subtle construction. Along the way we recover and discover numerous physical, combinatorial, polyhedral, and enumerative properties of Matryoshkas and cosmohedra, that we briefly outline in this introduction.

### 1.1 Combinatorial properties of the cosmohedron

The combinatorial structure of the cosmological wavefunction of \operatorname{Tr}(\phi^{3}) theory is described by the _Matryoshkas_ of an (n+2)-gon \pentagon_{n+2}, so these are the objects we wish to model geometrically in the cosmohedron. We call them Matryoshkas because their nesting structure is reminiscent of, though much richer than, that of the Matryoshka dolls created by S. Malyutin and V. Zvyozdochkin in Russia in the 1890s.

One way of constructing a Matryoshka, “from the inside out,” is to choose a subdivision of \pentagon_{n+2}, and then repeatedly choosing a group of adjacent outermost polygons and wrapping them into a larger polygon engulfing them. Alternatively, “from the outside in”, a Matryoshka is obtained recursively by choosing a subdivision of \pentagon_{n+2}, and then placing a Matryoshka inside each of its polygons.

Figure 1: An octagonal Matryoshka M and a maximal Matryoshka M^{\prime} containing it.

Naturally, we can think of a Matryoshka as a set of polygons. This gives a partial order on the set of Matryoshkas on \pentagon_{n+2} by _containment_: we have M\leq M^{\prime} if every polygon of M is also in M^{\prime}.

Arkani-Hamed, Figueiredo, and Vazão introduced an (n-1)-dimensional polytope called the _cosmohedron_ to encode this combinatorial structure; see Definition [4.15](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem15 "Definition 4.15. ‣ 4.2.4 Isomorphism with the Arkani–Hamed-Figueiredo-Vazão cosmohedron. ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron"). Our main result is a proof of the following statement from [[AHFV25](https://arxiv.org/html/2603.03425v3#bib.bibx8)].

###### Theorem 1.1.

The face lattice of the (n-1)-cosmohedron is anti-isomorphic to the poset of Matryoshkas of an (n+2)-gon.

### 1.2 Polyhedral properties of the cosmohedron

The search for a cosmohedron that models the combinatorics of Matryoshkas is analogous, and in fact closely related, to earlier constructions of polyhedra with a desired combinatorial structure. An important motivating example is the _associahedron_, whose faces are in bijection with the plane trees of a given size.

The cosmohedron is obtained from the associahedron through a careful chiseling procedure: at each vertex a_{T} corresponding to the binary tree T, we chisel a copy of the _bracket associahedron_\operatorname{Br}_{T}. The chiselings at the different vertices interact with each other, so they must be carried out very precisely. We must place the right realization of each \operatorname{Br}_{T} in the correct position, angle, and size, so that the different bracket associahedra glue correctly into one polytope with the desired combinatorics. The 3-cosmohedron is shown in Figure [2](https://arxiv.org/html/2603.03425v3#S1.F2 "Figure 2 ‣ 1.2 Polyhedral properties of the cosmohedron ‣ 1 Introduction ‣ Combinatorics of the cosmohedron"). The precise construction is given in Section [4](https://arxiv.org/html/2603.03425v3#S4 "4 The cosmohedron ‣ Combinatorics of the cosmohedron").

Figure 2: The 3-cosmohedron as a blowup of the 3-associahedron, and the Matryoshkas labeling some of its faces. 

Arkani–Hamed, Figueiredo, and Vazão [[AHFV25](https://arxiv.org/html/2603.03425v3#bib.bibx8)] proposed a construction of the cosmohedron \mathrm{AFV}_{n-1} as an (n-1)-dimensional polytope in \mathbb{R}^{(n+2)(n-1)/2}; see Definition [4.15](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem15 "Definition 4.15. ‣ 4.2.4 Isomorphism with the Arkani–Hamed-Figueiredo-Vazão cosmohedron. ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron"). This embedding gives the polytope an elegant inequality description, but it does not offer a direct path to prove its correctness, as we explain in Sections [1.4](https://arxiv.org/html/2603.03425v3#S1.SS4 "1.4 Elusive qualities of the cosmohedron ‣ 1 Introduction ‣ Combinatorics of the cosmohedron") and [6](https://arxiv.org/html/2603.03425v3#S6 "6 Elusive qualities of the cosmohedron ‣ Combinatorics of the cosmohedron").

Our approach is to construct a different embedding \operatorname{Cosmo}_{n-1} of the cosmohedron in (a hyperplane in) \mathbb{R}^{n} which is closer in spirit to the family of generalized permutahedra. To do so, we begin by building its normal fan \operatorname{CosmoFan}_{n-1} in \mathbb{R}^{n}/\mathbb{R}, and proving that it has the correct combinatorial structure of Matryoshkas; we do this in Definition [2.3](https://arxiv.org/html/2603.03425v3#S2.Thmtheorem3 "Definition 2.3. ‣ 2.2 The cosmohedral fan in terms of Matryoshkas ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron") and Theorem [3.6](https://arxiv.org/html/2603.03425v3#S3.Thmtheorem6 "Theorem 3.6. ‣ 3.3 The cosmohedral fan is a complete fan ‣ 3 Proof of correctness of the cosmohedral fan ‣ Combinatorics of the cosmohedron"). Next, we construct explicit coordinates for the vertices of \operatorname{Cosmo}_{n-1}, and prove that its normal fan is \operatorname{CosmoFan}_{n-1} in Theorem [4.11](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem11 "Theorem 4.11. ‣ 4.2.2 The cosmohedron: proof of correctness ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron"). Since we know the rays of the fan, we can readily use them to compute the inequality description of \operatorname{Cosmo}_{n-1} in Theorem [4.13](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem13 "Theorem 4.13. ‣ 4.2.3 Inequality description. ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron"). Finally, with inequality descriptions at hand for both of them, we find a linear isomorphism between \operatorname{Cosmo}_{n-1} in \mathbb{R}^{n} and \mathrm{AFV}_{n-1} in \mathbb{R}^{(n+2)(n-1)/2} in Theorem [4.16](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem16 "Theorem 4.16. ‣ 4.2.4 Isomorphism with the Arkani–Hamed-Figueiredo-Vazão cosmohedron. ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron"). This proves the combinatorial structure of the cosmohedron \mathrm{AFV}_{n-1} as conjectured in [[AHFV25](https://arxiv.org/html/2603.03425v3#bib.bibx8)].

### 1.3 Enumerative properties of the cosmohedron

In Section [5](https://arxiv.org/html/2603.03425v3#S5 "5 Face structure and enumeration ‣ Combinatorics of the cosmohedron") we turn to enumerative questions. The facets of the cosmohedron are in bijection with the non-trivial subdivisions of \pentagon_{n+2}, which are enumerated by the Hipparchus–Schröder numbers. This sequence is well understood, and is in fact one of the earliest combinatorial sequences, dating to the second century BC; for an interesting historical account, see [[Sta99](https://arxiv.org/html/2603.03425v3#bib.bibx37), [Sta97](https://arxiv.org/html/2603.03425v3#bib.bibx36)].

The vertices correspond to the maximal Matryoshkas; we prove that they are enumerated by sequence A177384 in the OEIS. We give a recursive formula for the sequence and a simple polynomial differential equation for the power series:

###### Theorem 1.2.

The generating function M(x)=x^{2}+2x^{3}+10x^{4}+72x^{5}+\cdots for the number of maximal Matryoshkas is characterized by the polynomial differential equation

M(x)=x^{2}+M(x)M^{\prime}(x).

Thus Matryoshkal numbers form one of the simplest sequences in the notoriously difficult family of _D-algebraic power series_[[Mel21](https://arxiv.org/html/2603.03425v3#bib.bibx30)]. For instance, we do not know how to prove the asymptotic behavior of this sequence; see Conjecture [5.7](https://arxiv.org/html/2603.03425v3#S5.Thmtheorem7 "Conjecture 5.7. ‣ 5.3 Asymptotic growth ‣ 5 Face structure and enumeration ‣ Combinatorics of the cosmohedron").

We also give an intriguing description for the f-vectors of cosmohedra: their generating function equals its own compositional inverse up to a simple transformation. Consider the f-polynomials of the cosmohedral fans f_{n}(t) as well as the simple transformation g_{n}(t)=((1+t)f_{n}(t)-1)/t. The coefficients of g_{n}(t) are obtained by adding consecutive coefficients of f_{n}(t):

\displaystyle f_{1}(t)=1,\,f_{2}(t)=1+2t,\,f_{3}(t)=1+10t+10t^{2},\,f_{4}(t)=1+44t+114t^{2}+72t^{3},\ldots
\displaystyle g_{1}(t)=1,\,g_{2}(t)=3+2t,\,g_{3}(t)=11+20t+10t^{2},\,g_{4}(t)=45+158t+186t^{2}+72t^{3},\ldots

###### Theorem 1.3.

The f-polynomials of cosmohedra form the unique sequence \{f_{n}(t)\}_{n\geq 1} such that the power series in (\mathbb{R}[t])[x]

\displaystyle x-f_{1}(t)x^{2}-f_{2}(t)x^{3}-f_{3}(t)x^{4}-\ldots
\displaystyle x+g_{1}(t)x^{2}+g_{2}(t)x^{3}+g_{3}(t)x^{3}+\ldots

are compositional inverses, where g_{n}(t)=((1+t)f_{n}(t)-1)/t for n\geq 1.

We also give a strikingly similar formula for the f-vector of the correlatron.

### 1.4 Elusive qualities of the cosmohedron

As mentioned above, Arkani–Hamed, Figueiredo, and Vazão’s cosmohedron \mathrm{AFV}_{n-1} is an (n-1)-dimensional polytope in \mathbb{R}^{(n+2)(n-1)/2}, given by a combinatorially rich system of inequalities; see Definition [4.15](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem15 "Definition 4.15. ‣ 4.2.4 Isomorphism with the Arkani–Hamed-Figueiredo-Vazão cosmohedron. ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron"). How might one prove its combinatorial structure? The standard methods don’t give a simple approach.

One might try to prove that \mathrm{AFV}_{n-1} has the correct incidences between vertices and facets. However, we did not have a formula for the vertices of \mathrm{AFV}_{n-1} before; and now that we do thanks to Proposition [4.8](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem8 "Proposition 4.8. ‣ 4.2.1 Vertex description ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron"), we see they are rather intricate. Proving the vertex-facet incidences would be a highly non-trivial matter.

To circumvent this issue, one might instead describe the normal fan of \mathrm{AFV}_{n-1}, and prove the corresponding wall-crossing inequalities, as described for example in [[ACEP20](https://arxiv.org/html/2603.03425v3#bib.bibx2), Lemma 2.10] or [[CLS24](https://arxiv.org/html/2603.03425v3#bib.bibx18), Theorems 6.1.5–6.1.7]. We know at the outset that this computation would be much more cumbersome than usual, because the normal fan is not simplicial. Furthermore, in this high dimensional embedding, the (n-1)-dimensional normal fan lives in a combinatorially intricate quotient of an \frac{1}{2}(n+2)(n-1)–dimensional space, so the description of the normal fan is far from unique. These obstacles make this computation very subtle as well.

These are the reasons why we proceeded by giving a different construction of the cosmohedron \operatorname{Cosmo}_{n-1} in \mathbb{R}^{n}, in the spirit of generalized permutahedra, and proved its combinatorial structure. Then, as we somewhat optimistically hoped at the outset, the polytopes \operatorname{Cosmo}_{n-1} and \mathrm{AFV}_{n-1} turned out to be linearly, and thus combinatorially equivalent.

However, proving the correctness of \operatorname{Cosmo}_{n-1} is still a delicate matter. Cosmohedra are related to several polytopes in the mathematical literature, some of which are known to be relevant in physical applications, and many of which are generalized permutahedra. These include permutahedra [[BS13](https://arxiv.org/html/2603.03425v3#bib.bibx14), [Sch13](https://arxiv.org/html/2603.03425v3#bib.bibx34)], associahedra [[CSZ15](https://arxiv.org/html/2603.03425v3#bib.bibx19), [Lod04](https://arxiv.org/html/2603.03425v3#bib.bibx29), [Sta63](https://arxiv.org/html/2603.03425v3#bib.bibx35), [Tam54](https://arxiv.org/html/2603.03425v3#bib.bibx38)], bracket associahedra [[BCL+22](https://arxiv.org/html/2603.03425v3#bib.bibx12), [LA22](https://arxiv.org/html/2603.03425v3#bib.bibx28)] (which are special cases of graph associahedra [[CD06](https://arxiv.org/html/2603.03425v3#bib.bibx15), [Dev09](https://arxiv.org/html/2603.03425v3#bib.bibx20), [Pos09](https://arxiv.org/html/2603.03425v3#bib.bibx32)]), permutoassociahedra [[CL23](https://arxiv.org/html/2603.03425v3#bib.bibx17), [Kap93](https://arxiv.org/html/2603.03425v3#bib.bibx25)], and permutonestohedra [[Gai15](https://arxiv.org/html/2603.03425v3#bib.bibx24)], among others. There are numerous methods to prove the combinatorial structure of such polytopes.

In Section [6](https://arxiv.org/html/2603.03425v3#S6 "6 Elusive qualities of the cosmohedron ‣ Combinatorics of the cosmohedron") we explain that the cosmohedron \operatorname{Cosmo}_{n-1} is more subtle than its predecessors in several ways: it is not very symmetric, it is not simple, it depends very delicately on the chiseling, we do not know how to express it as a deformation of an easy polytope, we do not know how to express it as a Minkowski sum of easy polytopes, and its faces factor combinatorially but not geometrically. For these reasons, proving its combinatorial structure requires new ideas, as described above.

### 1.5 Cosmological properties of the cosmohedron

Let us close this introduction by briefly explaining the physical origins of this project. The content of this section is not essential to the mathematical understanding of this paper, but it is essential to the discovery of these mathematical objects, and it gives rise to many interesting new directions.

The past decade has seen the discovery of startling and deep new connections between fundamental physics and new mathematical structures in combinatorics and geometry. The story began with the most basic and fundamental physical processes in nature–the collision and scattering of elementary particles. Because of quantum mechanics, the outcome of any given scattering process is uncertain; only the probability for different states can be predicted. Given an initial state i and a final state f, the probability is |{\cal A}_{i,f}|^{2} where {\cal A}_{i,f} is the _scattering amplitude_. The textbook method for computing these amplitudes is to sum over _Feynman diagrams_ representing the different ways in which the collision events can take place in space and time.

The connection to combinatorics and geometry can most easily and vividly be illustrated in an especially useful toy model, the so-called _Tr(\phi^{3})_ theory, for “colored” scalar particles represented by N\times N square matrices. In the simplest limit (“tree-level”), the amplitude turns into the sum over all planar trivalent (known in physics as _cubic_) tree diagrams, or what is the same, a sum over the triangulations of an n-gon.

A priori to a physicist, there is no obvious combinatorics here: we simply have a collection of diagrams that have to be collected and summed over, in no particular order, to arrive at the final result. But a combinatorialist knows that there is something very special about this particular collection of objects: the triangulations of an n-gon are captured by the vertices of a polytope, the associahedron.

This link turns out to be the tip of a large iceberg of much larger significance to physics beyond “merely” the combinatorics of diagrams. A particular embedding of the associahedron allows us to extract mathematical expressions for the scattering amplitude from the geometry[[AHBHY18a](https://arxiv.org/html/2603.03425v3#bib.bibx4)]. This particular embedding makes it natural to think of the associahedron as a Minkowski sum of the simplest possible pieces – a collection of simplices – and quite remarkably, this fact is intimately related to (and indeed lets us naturally discover from the bottom-up) the generalization from particle scattering to scattering processes in string theory. Over the past few years, these developments have led to the discovery of a number of surprising connections between seemingly totally different physical theories[[AHCD+24](https://arxiv.org/html/2603.03425v3#bib.bibx6)] and, almost as an incidental byproduct, have provided new ways to simplify enormously complicated Feynman diagram calculations, sometimes into a few lines.

The connection to physics also led to the discovery of new mathematical objects. For instance, the associahedron captures all tree diagrams, but what about diagrams that can be drawn on surfaces with arbitrary topology? The consistency of the physical picture demands that geometric objects encoding such diagrams must exist, and they have indeed been described in the physics literature [[AHFS+25](https://arxiv.org/html/2603.03425v3#bib.bibx7)].

This emerging connection between physics and combinatorics suggests a radical new picture for the laws of physics, where the notion of particle evolution in spacetime is not primary. Instead, the unifying combinatorial-geometric objects that control how all spacetime processes are “glued together” play a more fundamental role. This dovetails nicely with an expectation widely shared by most theoretical physicists that the notion of spacetime cannot be fundamental but must eventually be supplanted by more elementary and abstract ideas. The new connection to combinatorics suggests a concrete strategy for pursuing this vision.

Nowhere is the need for a new picture of physics supplanting spacetime more urgent than in cosmology. Amongst other things, the notion of time itself breaks down at the big bang singularity, so the standard picture of time evolution must likely be abandoned.

Cosmology also offers us a rich set of physical observables to compute, and in principle measure, that turn out to be interesting generalizations of scattering amplitudes. It is widely believed that the very early history of the universe involved a period of exponentially rapid expansion known as _inflation_, and that quantum fluctuations in this early inflationary phase were blown up to enormous scales, even larger than our observable universe, giving rise to all the “clumpiness” in matter and all the interesting structures we see in the universe today. The details of these quantum fluctuations are reflected in the statistics of the distribution of mass and energy in the universe today and are known as _cosmological correlation functions_. The cosmological correlators produced during inflation are computed in much the same way as particle scattering amplitudes, drawing diagrams for all possible ways the correlations could have been produced. But there is an important difference: while for scattering processes a single diagram represents a single elementary process, in cosmology, even a single diagram has much more structure, precisely having to do with the different time-ordering of “histories” a given diagram can represent.

If we return to the toy Tr(\phi^{3}) theory, to compute a cosmological correlator, a given diagram must be further decorated by all the possible _bracketings_ of the graph, or equivalently, what we will call the possible _Matryoshkas_ associated with the corresponding triangulation. This arises in a very natural way when solving the Schrödinger equation H\Psi=E\Psi to compute the wavefunction \Psi(\phi) for \phi, from which we determine probabilities and correlation functions. The Schrödinger equation for \Psi(\phi)={\rm exp}[\-\psi(\phi)] can be solved in perturbation theory. The exponential relating \Psi(\phi) to \psi(\phi) reflects the familiar exponential relationship between all graphs and connected graphs. Every connected cubic diagram with n external legs, where the external legs carry spatial momentum \vec{k}_{i}, makes a contribution to the coefficient of \phi(\vec{k}_{1})\cdots\phi(\vec{k}_{n}) in the Taylor expansion of \psi(\phi) in powers of \phi. The contribution of a diagram to \psi depends on “energies” entering each subgraph – this is the total sum of the magnitude of the momenta entering the subgraph. In Figure [3](https://arxiv.org/html/2603.03425v3#S1.F3 "Figure 3 ‣ 1.5 Cosmological properties of the cosmohedron ‣ 1 Introduction ‣ Combinatorics of the cosmohedron"), these energies are given by perimeters of subpolygons.

![Image 1: Refer to caption](https://arxiv.org/html/2603.03425v3/Asset_4.png)

Figure 3: The contribution of a tree graph (or equivalently a triangulation) to the wavefunction is determined by the bracketings of the graph (or equivalently the Matryoshkas of the triangulation). Here k_{i,j} denotes the length between vertices i and j, which corresponds to the length of the momenta flowing through the edges of the graph, and \mathcal{P} stands for the perimeter of the subpolygons entering the Matryoshka.

Written in terms of \psi, the Schrödinger equation turns into a simple recursion relation: the contribution from a given graph can be computed by first writing a factor of the inverse of the sum over all energies on the vertices, and then summing over all the ways of deleting an edge from the graph, and adding the energy of the edge to the vertices it connects to in the smaller graph with the edge deleted. We can run this recursion until there are no more edges to be cut, and every term encountered corresponds to a bracketing of the graph – where to each bracket subgraph we associate the factor given by the total energy entering that subgraph. Phrased in the language of the triangulations, each bracket corresponds to a sub-polygon in the triangulation, to which we associate a factor given by its inverse perimeter. The terms contributing to the wavefunction are then associated with a maximal collection of nested sub-polygons in the triangulation, corresponding to a Matryoshka; see Figure [3](https://arxiv.org/html/2603.03425v3#S1.F3 "Figure 3 ‣ 1.5 Cosmological properties of the cosmohedron ‣ 1 Introduction ‣ Combinatorics of the cosmohedron").

Given the close connection between scattering amplitudes and cosmological correlators, physics suggested that some analog of the associahedron must exist, capturing not the combinatorics of all triangulations, but the factorially larger combinatorics of all possible Matryoshkas. This was the motivation leading to the discovery and definition of cosmohedra, which are polytopes precisely realizing this picture.

Our main goal in this paper is to better understand the combinatorics and geometry of cosmohedra, provide proofs of many of their claimed properties, and also describe them in a way that makes more natural contact with similar objects that have recently been studied in combinatorics.

### 1.6 Outlook: \mathcal{X} in Y polytopes and related constructions in physics

Apart from its intrinsic mathematical interest, our treatment makes it clear that cosmohedra belong to a wider class of objects whose relevance to physics should extend far beyond cosmology. We will see that cosmohedra generalize to what we can call _\mathcal{X} in Y_ polytopes, where we begin with a polytope Y and “chisel” it in a way controlled by a family \mathcal{X} of polytopes to get a new _chiseled polytope_. At the level of fans, the normal fan of a (simple) polytope Y is canonically refined using the normal fans of the polytopes in \mathcal{X}. We expect these \mathcal{X} in Y polytopes to be relevant to many settings in physics, where, in addition to summing over all diagrams, we must keep track of various properties of subgraphs in each diagram.

An especially important example of this arises in the context of scattering amplitudes beyond “tree” level. In diagrams with loops, the internal momenta are not determined uniquely and must be integrated over. This is associated with a rich pattern of possible divergences when loop momenta become large (“ultraviolet”) or small (“infrared”). For an individual Feynman diagram, these divergences are controlled by the Newton polytopes of certain graph polynomials – the so-called _Symanzik polynomials_ –associated with the graph. In [[AHCD+24](https://arxiv.org/html/2603.03425v3#bib.bibx6)], all Feynman diagrams are naturally realized as cones glued together into a complete “Feynman fan”, but this raises a question posed in [[AHCD+24](https://arxiv.org/html/2603.03425v3#bib.bibx6)]: can the combinatorics of the “Symanzik polytopes” associated with any graph be encoded in some way inside the full Feynman fan that glues all diagrams together? This is exactly what \mathcal{X} in Y polytopes allows us to do. In the Appendix, we sketch this construction in the first interesting case: the Symanzik \mathcal{U}-polynomials at one-loop level in the Tr(\phi^{3}) theory.

## 2 Matryoshkas, bracketed trees, and the cosmohedral fan

### 2.1 Matryoshkas

Consider a convex (n+2)-gon P. A _subdivision_\mathcal{S} is a set of sub-polygons whose union is P and whose interiors are pairwise disjoint. A subdivision \mathcal{S} can be naturally identified with the corresponding set \mathcal{C} of non-crossing diagonals of P. We say \mathcal{S}_refines_\mathcal{S}^{\prime}, and write \mathcal{S}\geq\mathcal{S}^{\prime}, if the following three equivalent conditions hold: every polygon in \mathcal{S} is contained in a polygon in \mathcal{S}^{\prime}, every polygon in \mathcal{S}^{\prime} is a union of polygons in \mathcal{S}, and the corresponding sets of chords satisfy \mathcal{C}\supseteq\mathcal{C}^{\prime}.

A _Matryoshka_ or _Russian doll_ M is a set of sub-polygons of P such that for any non-minimal polygon X in M, the maximal subpolygons (M_{<X})_{\text{max}} of X that are in the Matryoshka M subdivide X. In particular, the set M_{\text{min}} of minimal polygons in M form the _underlying subdivision_ of M. We say M _refines_ M^{\prime}, and write M\geq M^{\prime}, if the sets of polygons satisfy M\supseteq M^{\prime}.

For example, Figure [4](https://arxiv.org/html/2603.03425v3#S2.F4 "Figure 4 ‣ 2.1 Matryoshkas ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron") shows two maximal Matryoshkas M_{1} and M_{2} on a hexagon. Removing a quadrilateral from M_{1} gives a Matryoshka smaller than it. Removing a single triangle from M_{1} does not give a non-Matryoshka, because the corresponding quadrilateral is not subdivided by its maximal subpolygons.

Figure 4: Two maximal Matryoshkas M_{1} and M_{2} on a hexagon. 

Recursively, a Matryoshka is either:

\bullet M=\{P\}, the full polygon, or

\bullet M=\mathcal{S}\cup\bigcup_{S\in\mathcal{S}}M_{S}: a subdivision \mathcal{S} of P and a Matryoshka M_{S} on each part S of \mathcal{S}.

###### Definition 2.1.

The _poset of Matryoshkas_\operatorname{Mat}_{n-1} is the partially ordered set of Matryoshkas on the (n+2)-gon ordered by containment, with an additional maximum element \widehat{1}.

### 2.2 The cosmohedral fan in terms of Matryoshkas

Let us give the vertices of our (n+2)-gon P=\pentagon_{n+2} the labels 0,1,\ldots,n+1 counterclockwise, with the edge e joining 0 and n+1 placed horizontally at the top.

Consider a triangulation (or subdivision) \mathcal{S} of P.

Triangle (or polygon) labels. We have a canonical labeling of the n triangles of \mathcal{S} with the vertex labels 1,\ldots,n: We give the label v to the first triangle of \mathcal{S} that one encounters on a straight path from vertex v to edge e. (In a subdivision, a k-gon gets k-2 such labels.)

Diagonal labels. Each diagonal of \mathcal{S} separates two triangles labeled v and w (or polygons labeled V and W), with triangle v between triangle w and edge e; we give that diagonal the label vw^{\perp} (or VW^{\perp}).

Now consider a Matryoshka M of P, with underlying triangulation (or subdivision) \mathcal{S}.

Polygon labels. Each diagonal vw^{\perp} is contained inside a unique smallest polygon of M, which we denote \pentagon_{M}(vw^{\perp}). When M is maximal, this is a bijection between the diagonals of \mathcal{S} and the non-minimal polygons of M.

Containment poset. The diagonals vw^{\perp} are the elements of the _containment poset_\tau of M, whose order relations are given by containment of the corresponding polygons: ij^{\perp}\geq kl^{\perp} if \pentagon_{M}(ij^{\perp})\supseteq\pentagon_{M}(kl^{\perp}). This poset is a tree, which is binary if M is maximal.

Matryoshkal cone. To the Matryoshka M we assign the cone

\operatorname{cone}(M)=\{x\in\mathbb{R}^{n}\,:\,x_{i}-x_{j}\geq x_{k}-x_{l}\geq 0\text{ when $\pentagon_{M}(ij^{\perp})\supseteq\pentagon_{M}(kl^{\perp})$}\}(1)

naturally associated to the containment poset.

Figure 5: The Matryoshka M_{1} and its containment poset \tau_{1}.

###### Example 2.2.

In the Matryoshka M_{1} of Figure [5](https://arxiv.org/html/2603.03425v3#S2.F5 "Figure 5 ‣ 2.2 The cosmohedral fan in terms of Matryoshkas ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron") we have \pentagon(21^{\perp})\subset\pentagon(32^{\perp}) and \pentagon(34^{\perp})\subset\pentagon(32^{\perp}), so \operatorname{cone}(M)=\{x\in\mathbb{R}^{4}\,:\,x_{3}-x_{2}\geq\{x_{2}-x_{1},x_{3}-x_{4}\}\geq 0\}.

###### Definition 2.3.

The _cosmohedral fan_ is

\operatorname{CosmoFan}_{n-1}=\{\operatorname{cone}(M)\,:\,M\text{ is a Matryoshka of the $(n+2)$-gon}\}.

We will prove in Section [3](https://arxiv.org/html/2603.03425v3#S3 "3 Proof of correctness of the cosmohedral fan ‣ Combinatorics of the cosmohedron") that this is indeed a fan, and that it has the desired combinatorial structure:

###### Theorem 2.4.

The cosmohedral fan is a complete fan in \mathbb{R}^{n}. The map M\mapsto\operatorname{cone}(M) is an order-preserving bijection between the poset \operatorname{Mat}_{n-1} of Matryoshkas of an (n+2)-gon and the face lattice L(\operatorname{CosmoFan}_{n-1}) of the cosmohedral fan.

### 2.3 Bracketed trees and the cosmohedral fan

Let us give a different combinatorial description of Matryoshkas, which gives rise to a description of the cosmohedral fan in terms of _bracketed trees_.

#### 2.3.1 Preliminaries: plane trees, triangulations, and the associahedral fan.

A _plane tree_ is a rooted tree where the children of each node are linearly ordered from left to right. The plane trees on n+1 leaves form a poset \operatorname{Tree}_{n-1} where T\leq T^{\prime} if T can be obtained from T^{\prime} by contracting edges. We include an additional maximum element \widehat{1} that makes this poset into a lattice.

The coatoms of this poset are the _binary trees_, which have n internal nodes and n+1 leaves. It is sometimes convenient to erase the leaves of a binary tree, obtaining a plane tree with n vertices.

If \mathcal{S} is a subdivision of \pentagon_{n+2}, then the _dual graph_ T(\mathcal{S}) has a vertex for each polygon of \mathcal{S} and an edge between each pair of polygons that share an edge. This is a tree rooted at the polygon adjacent to the top edge e. Figure [9](https://arxiv.org/html/2603.03425v3#S2.F9 "Figure 9 ‣ 2.3.3 The cosmohedral fan in terms of bracketed trees. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron") shows the binary trees dual to the subdivisions of Matryoshkas M_{1} and M_{2}.

If \mathcal{S} is a triangulation then T(\mathcal{S}) is a binary tree, and this is a bijection between the triangulations of \pentagon_{n+2} and the trees in \operatorname{Tree}_{n-1}; the number of such objects is the Catalan number C_{n}=\frac{1}{n+1}\binom{2n}{n}.

The internal vertices of a plane tree T with n+1 leaves are naturally labeled with sets that form a partition of [n]: we give the gaps between the leaves the labels 1,2,\ldots,n from left to right, and label each internal node of T with the gap(s) that are visible from it. We regard T as the Hasse diagram of a pre-poset on [n] whose root is the maximum element, and define

\operatorname{cone}(T)=\{x\in\mathbb{R}^{n}/\mathbb{R}\mathbf{1}\,:\,x_{i}\geq x_{j}\text{ for each }i\geq j\text{ in }T_{\operatorname{int}}\}

where \mathbf{1}=(1,\ldots,1). In particular x_{i}=x_{j} in \operatorname{cone}(T) whenever i and j label the same vertex of T.

###### Definition 2.5.

The _associahedral fan_ is

\operatorname{AssocFan}_{n-1}=\{\operatorname{cone}(T)\,:\,T\text{ is a plane tree on $n+1$ leaves}\}.

Figure [6](https://arxiv.org/html/2603.03425v3#S2.F6 "Figure 6 ‣ 2.3.1 Preliminaries: plane trees, triangulations, and the associahedral fan. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron") shows the associahedral fan \operatorname{AssocFan}_{3}.

![Image 2: Refer to caption](https://arxiv.org/html/2603.03425v3/Fan_Labels2.png)

Figure 6: The three-dimensional associahedral fan.

###### Theorem 2.6.

[[Lod04](https://arxiv.org/html/2603.03425v3#bib.bibx29)] The associahedral fan \operatorname{AssocFan}_{n-1} is a complete fan in \mathbb{R}^{n}. The map T\mapsto\operatorname{cone}(T) is an order-preserving bijection between the poset \operatorname{Tree}_{n-1} of plane trees with n+1 leaves and the face lattice L(\operatorname{AssocFan}_{n-1}) of the associahedral fan.

#### 2.3.2 Preliminaries: graph bracketings and bracket associahedral fans.

Let H=(V,E) be a connected graph. A _bracket_ of H is a connected subgraph of H with at least one edge. We identify it with its set of edges. A _bracketing_ of H is a collection B of brackets such that any two brackets are either nested (one is contained in the other) or disjoint (they have no edges or vertices in common). We require every bracketing to contain the full graph H.

The _poset of bracketings_\operatorname{Br}(H) is the set of bracketings of H ordered by containment, with an additional maximum element \widehat{1}.

In a bracketing B, each edge e is contained in a unique smallest bracket of B, which we denote \operatorname{bracket}_{B}(e) or simply \operatorname{bracket}(e). This defines a pre-poset \operatorname{tree}(B) on the edge set E, where e\geq_{B}f if e is _outermore_ than f in the sense that \operatorname{bracket}(e)\supseteq\operatorname{bracket}(f); or equivalently, every bracket containing e also contains f. Because every bracket is contained in a unique smallest bracket, the poset \operatorname{tree}(B) is a tree rooted at a unique maximal element. We obtain a cone

\operatorname{cone}_{H}(B)=\{y\in\mathbb{R}^{E(H)}\,:\,y_{e}\geq y_{f}\text{ for each }e\geq_{B}f\}.

Its dimension is |B|.

###### Definition 2.7.

The _bracket associahedral fan_ of a graph H is

\operatorname{BrAssocFan}(H)=\{\operatorname{cone}_{H}(B)\,:\,B\text{ is a bracketing of $H$}\}

Figure [7](https://arxiv.org/html/2603.03425v3#S2.F7 "Figure 7 ‣ 2.3.2 Preliminaries: graph bracketings and bracket associahedral fans. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron") illustrates the bracket associahedral fan of a path of length 4 with edges e,f,g. The fan is three-dimensional, with all cones containing the line y_{e}=y_{f}=y_{g}, so we display a 2-dimensional slice.

###### Theorem 2.8.

[[AHFV25](https://arxiv.org/html/2603.03425v3#bib.bibx8), [CD06](https://arxiv.org/html/2603.03425v3#bib.bibx15)] The bracket associahedral fan \operatorname{BrAssocFan}(H) is a complete fan in \mathbb{R}^{E(H)}. The map B\mapsto\operatorname{cone}_{H}(B) is an order-preserving bijection between the poset of bracketings \operatorname{Br}(H) and the face lattice L(\operatorname{BrAssocFan}(H)) of the bracket associahedral fan.

Figure 7: The bracket associahedral fan of a path.

Bracket associahedral fans and graph associahedral fans. The _line graph_ L(H) of H=(V,E) is the graph of edge adjacencies of H: it has a vertex for each edge e of H, and the vertices corresponding to edges e,f are connected in L(H) if e and f share a vertex in H. The bracketings of H are in order-preserving bijection with the tubings of L(H) in the sense of [[CD06](https://arxiv.org/html/2603.03425v3#bib.bibx15)], so

\{\text{bracket associahedral fans}\}\subsetneq\{\text{graph associahedral fans}\}.

Figure 8: Bracketings of H correspond to tubings of L(H).

The graph associahedral fan of the claw K_{1,3} – studied and depicted in [[ARW06](https://arxiv.org/html/2603.03425v3#bib.bibx11), Figure 3] in its manifestation as the Dynkin diagram D_{4} – is a three-dimensional fan with 10 rays, corresponding to the 10 tubes of the graph. The reader may verify that this is the only graph on four vertices that has 10 tubes. Since the claw is not a line graph, this is not a bracket associahedral fan. It is then natural to ask:

Which graph associahedral fans are bracket associahedral fans?

#### 2.3.3 The cosmohedral fan in terms of bracketed trees.

A _bracketed tree_(T,B) consists of a plane tree T and a bracketing B of T. Let us define the cone of a bracketed tree (T,B) with n+1 leaves to be

\operatorname{cone}(T,B)=\{x\in\mathbb{R}^{n}\,:\,x_{i}-x_{j}\geq x_{k}-x_{l}\geq 0\text{ when $\operatorname{bracket}_{B}(ij)\supseteq\operatorname{bracket}_{B}(kl)$}\}.

Note that the inequalities of \operatorname{cone}(T,B) correspond to the order relations of the tree poset \operatorname{tree}(B), as illustrated in Figure [9](https://arxiv.org/html/2603.03425v3#S2.F9 "Figure 9 ‣ 2.3.3 The cosmohedral fan in terms of bracketed trees. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron").

The goal of this section is to show that the cones of the bracketed trees on n+1 leaves are precisely the cones of the cosmohedral fan \operatorname{CosmoFan}_{n-1}. We begin with a technical lemma.

###### Lemma 2.9.

For any bracketed tree (T,B) we have

\operatorname{cone}(T,B)\subseteq\operatorname{cone}(T),\qquad\dim\operatorname{cone}(T,B)=\dim\operatorname{cone}(T)=|V(T)_{int}|+1.

###### Proof.

The inclusion and the formula for \dim\operatorname{cone}(T) are clear from the definitions. To prove that \operatorname{cone}(T,B) has the same dimension, we first claim that for x\in\operatorname{cone}(T,B) we have that x_{i}=x_{j} if i and j label the same vertex of T. To see this, note that for an offspring k we have \operatorname{bracket}_{B}(ik)=\operatorname{bracket}_{B}(jk) so x_{i}-x_{k}=x_{j}-x_{k} (or analogously for the parent k, if they do not have offspring). Setting these equalities aside, the remaining inequalities of \operatorname{cone}(T,B) are given by the acyclic poset \operatorname{tree}(B), so they imply no further equalities. ∎

Figure 9: Two Matryoshkas, their corresponding bracketed trees, and their cones.

The bijection with Matryoshkas. For a Matryoshka M, let \mathcal{S}=M_{\min} be the set of inclusion-minimal polygons in M, which form the underlying subdivision of M.

1.   1.
The _rooted planar tree_ T_{M} of M is the dual graph of \mathcal{S}, rooted at the top polygon containing edge e.

2.   2.
Each non-minimal polygon P in M is a disjoint union of minimal polygons in M_{\min}, whose dual graph is a bracket of the tree T_{M}, denoted \operatorname{bracket}(P). The _bracketing_ B_{M} of M is B_{M}=\{\operatorname{bracket}(P)\,:\,P\in M_{not\,min}\}.

###### Lemma 2.10.

The map M\mapsto(T_{M},B_{M}) is a bijection between the Matryoshkas on an (n+2)-gon and the bracketed trees \operatorname{BrTree}_{n-1} with n+1 leaves.

###### Proof.

The Matryoshka condition on M implies that B_{M} is indeed a bracketing of T_{M}. Conversely, for any pair (T,B), the plane tree T is dual to a subdivision \mathcal{S} of the (n+2)-gon, and each bracket in B is dual to a polygon P that is the union of polygons in \mathcal{S}. The bracketing condition on B implies that these polygons form a Matryoshka. ∎

###### Lemma 2.11.

Under the bijection of Lemma [2.10](https://arxiv.org/html/2603.03425v3#S2.Thmtheorem10 "Lemma 2.10. ‣ 2.3.3 The cosmohedral fan in terms of bracketed trees. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron"), the cone of the Matryoshka M equals the cone of the bracketed tree (T_{M},B_{M}):

\operatorname{cone}(M)=\operatorname{cone}(T_{M},B_{M}).

###### Proof.

This follows readily from the definitions. ∎

###### Corollary 2.12.

We have

\operatorname{CosmoFan}_{n-1}=\{\operatorname{cone}(T,B)\,:\,(T,B)\text{ is a bracketed tree with $n+1$ leaves}\}.

The poset of bracketed trees. Consider a graph G and a bracketing B. If \beta is a bracket of B, the _deletion_ is the bracketing B-\beta is a bracketing of G. If \mu is a minimal bracket of B, the _contraction_ is the bracketing B/\mu:=\{\beta-\mu\,:\,\beta\in B,\beta\neq\mu\} of the tree T/\mu where edge \mu is contracted – so the edge is deleted and its vertices are identified; see Figure [10](https://arxiv.org/html/2603.03425v3#S2.F10 "Figure 10 ‣ 2.3.3 The cosmohedral fan in terms of bracketed trees. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron"). If E is a downset in (the tree pre-poset of) B, we can define B/E by successively contracting the brackets in E in a non-decreasing order.

Figure 10: A contraction of a bracketing.

The _poset of bracketed plane trees_ with n+1 leaves, denoted \operatorname{BrTree}_{n-1}, is defined by deletion and contraction as follows. Locally, we have a cover relation (T,B)\lessdot(T^{\prime},B^{\prime}) if either:   
(i) T=T^{\prime} and B=B^{\prime}-\beta^{\prime} for a bracket \beta^{\prime} of B^{\prime}, or   
(ii) T=T^{\prime}/\mu^{\prime} and B=B^{\prime}/\mu^{\prime} for a minimal bracket \mu^{\prime} of B^{\prime}.   
Globally, we have a poset relation (T,B)\leq(T^{\prime},B^{\prime}) if T=T^{\prime}/E for a set of edges E that is a downset of the tree poset of B^{\prime}, and B\subseteq B^{\prime}/E. We augment this poset with a maximum element \widehat{1}.

###### Proposition 2.13.

The bijection M\mapsto(T_{M},B_{M}) of Lemma [2.10](https://arxiv.org/html/2603.03425v3#S2.Thmtheorem10 "Lemma 2.10. ‣ 2.3.3 The cosmohedral fan in terms of bracketed trees. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron"), is a poset isomorphism between the lattice of Matryoshkas \operatorname{Mat}_{n+2} and the poset of bracketed trees \operatorname{BrTree}_{n-1}.

###### Proof.

In a cover relation M\lessdot M^{\prime} of Matryoshkas, M is obtained from M^{\prime} by either:   
(i) removing a non-minimal polygon P, which corresponds to leaving T_{M}=T_{M^{\prime}} fixed and making B_{M}=B_{M^{\prime}}-\operatorname{bracket}(P), or   
z (ii) removing the minimal polygons P_{1},\ldots,P_{k} that subdivided a next-to-minimal polygon P in M^{\prime}, so P becomes minimal in M; this corresponds to contracting the minimal \operatorname{bracket}(P) to get (T_{M},B_{M})=(T_{M^{\prime}}/\operatorname{bracket}(P),B_{M^{\prime}}/\operatorname{bracket}(P)).

This proves the poset isomorphism. ∎

## 3 Proof of correctness of the cosmohedral fan

### 3.1 Positive bracket associahedral fans

The _positive bracket associahedral fan_ of a graph H is the intersection

\operatorname{BrAssocFan}_{\geq 0}(H)=\operatorname{BrAssocFan}(H)\cap\mathbb{R}_{\geq 0}^{E(H)}

of the bracket associahedral fan with the positive orthant, with the inherited polyhedral structure.

Figure 11: The positive bracket associahedral fan of a path.

###### Example 3.1.

Figure [11](https://arxiv.org/html/2603.03425v3#S3.F11 "Figure 11 ‣ 3.1 Positive bracket associahedral fans ‣ 3 Proof of correctness of the cosmohedral fan ‣ Combinatorics of the cosmohedron") illustrates the positive bracket associahedral fan of a path of length 4 with edges e,f,g; it is the intersection of Figure [7](https://arxiv.org/html/2603.03425v3#S2.F7 "Figure 7 ‣ 2.3.2 Preliminaries: graph bracketings and bracket associahedral fans. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron") – which we should recall is three-dimensional – with the positive orthant. The maximal cones are in bijection with those of Figure [7](https://arxiv.org/html/2603.03425v3#S2.F7 "Figure 7 ‣ 2.3.2 Preliminaries: graph bracketings and bracket associahedral fans. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron"), but there are new lower-dimensional cones, corresponding to the bracketings B of the contractions G of the path. We label a few of the faces with a picture of the _pre-bracketing_(G,B), and with the list B of brackets, which is enough to recover G. We now prove this description holds in general.

A _pre-bracketing_ is a bracketing of a contraction of H; that is, a pair (G,B) consisting of a contraction G of H and a bracketing B of G. When H is fixed, we simply write B instead of (G,B), since G is determined automatically as the maximum bracket in B.

The _poset of pre-bracketings_\operatorname{PreBr}(H) is again defined by deletion and contraction: B\lessdot B^{\prime} if B is obtained from B^{\prime} by either (i) removing a bracket or (ii) contracting a minimal bracket.

We define the _positive cone of a pre-bracketing_ B of H to be

\operatorname{cone}_{H}(B)_{\geq 0}=\{y\in\mathbb{R}^{E(H)}\,:\,y_{e}\geq y_{f}\geq 0\text{ for each }e\geq_{B}f\,,\,\,\,y_{e}=0\text{ for each }e\text{ contracted}\}.

###### Proposition 3.2.

For any graph H, the map B\mapsto\operatorname{cone}_{H}(B)_{\geq 0} is an order-preserving bijection between the poset \operatorname{PreBr}(H) of pre-bracketings of H and the face lattice   
L(\operatorname{BrAssocFan}_{\geq 0}(H)) of the positive bracket associahedral fan of H.

###### Proof.

A face of the positive bracket associahedral fan is obtained by intersecting a face \operatorname{cone}_{H}(B) of \operatorname{BrAssocFan}(H) for a bracketing B with a face \mathbb{R}_{\geq 0}^{E} of \mathbb{R}_{\geq 0}^{E(H)} for a subset E\subseteq E(H). Imposing the equations y_{f}=0 for f\notin E on \operatorname{cone}_{H}(B) also implies y_{f^{\prime}}=0 for any f^{\prime}<_{B}f; so if we let \overline{F} be the downset of (the tree poset of) B generated by F=E(H)-E, then

\operatorname{cone}_{H}(B)\cap\mathbb{R}_{\geq 0}^{E}=\operatorname{cone}_{H}(B/\overline{F})_{\geq 0}.(2)

Every pre-bracketing of H arises in this way, proving the description of the cones of \operatorname{BrAssoc}_{\geq 0}(H).

Now, for a pre-bracketing B of H, a facet of \operatorname{cone}_{H}(B)_{\geq 0} is obtained by either (i) setting y_{e}=y_{f} for e\gtrdot f in \operatorname{tree}(B), which has the effect of deleting the bracket of f, or (ii) setting y_{f}=0 for a minimal f in \operatorname{tree}(B), which has the effect of contracting the bracket of f. The poset isomorphism follows. ∎

### 3.2 The cosmohedral fan and positive bracket associahedral fans

Now we describe the cosmohedral fan in terms of positive bracket associahedral fans.

###### Lemma 3.3.

Let T be a plane rooted tree with n+1 leaves, whose internal vertices are canonically labelled by [n]. Label each edge in the order ij where i is closer to the root than j. We have inverse maps:

\displaystyle f_{T}:\mathbb{R}^{E(T)}\displaystyle\overset{\cong}{\longrightarrow}\mathbb{R}^{n}/\mathbb{R}\mathbf{1}\displaystyle\qquad f_{T}^{-1}:\qquad\qquad\mathbb{R}^{n}/\mathbb{R}\mathbf{1}\displaystyle\overset{\cong}{\longrightarrow}\mathbb{R}^{E(T)}
\displaystyle e_{ij}\displaystyle\longmapsto\sum_{k\in(T-ij)_{i}}e_{k}+\mathbb{R}\displaystyle\qquad(x_{v})_{v\in V(T)}+\mathbb{R}\mathbf{1}\displaystyle\longmapsto(x_{i}-x_{j})_{ij\in E(T)}.

where (T-ij)_{i} and (T-ij)_{j} are the connected components of the forest T-ij containing i and j respectively. Furthermore,

1.   1.
the map f_{T}:\mathbb{R}^{E(T)}\overset{\cong}{\longrightarrow}\mathbb{R}^{n}/\mathbb{R}\mathbf{1} is an isomorphism of vector spaces,

2.   2.
the restriction f_{T}:\mathbb{R}_{\geq 0}^{E(T)}\overset{\cong}{\longrightarrow}\operatorname{cone}(T) is an isomorphism of cones, and

3.   3.
the restriction f_{T}:\operatorname{cone}_{T}(B)_{\geq 0}\overset{\cong}{\longrightarrow}\operatorname{cone}(T,B) is an isomorphism of cones for any bracketing B of tree T.

4.   4.
the restriction f_{T}:\operatorname{cone}_{T}(B)_{\geq 0}\overset{\cong}{\longrightarrow}\operatorname{cone}(T(B),B) is an isomorphism of cones for any pre-bracketing B of tree T, where T(B)\leq T is the tree supporting B.

###### Proof.

First, we check that these maps are indeed inverse to each other. For each basis vector e_{ij} corresponding to edge ij of T, the vector f_{T}(e_{ij}) has coordinates equal to 1 on component (T-ij)_{i} and 0 on component (T-ij)_{j}. Therefore f_{T}^{-1}(f_{T}(e_{ij})) is only non-zero on the edge ij – whose vertices lie in different components – where it equals 1. Thus f_{T}^{-1}(f_{T}(e_{ij}))=e_{ij} indeed. Since both spaces are (n-1)-dimensional, this means f_{T} is an isomorphism and f_{T}^{-1} is its inverse, proving 1.

To show 2, notice that y=f_{T}^{-1}(x) is nonnegative if and only if x satisfies the inequalities of \operatorname{cone}(T). For 3., observe that y satisfies the additional inequalities y_{ij}\geq y_{kl}\geq 0 of \operatorname{cone}_{T}(B)_{\geq 0} if and only if x satisfies the inequalities x_{i}-x_{j}\geq x_{k}-x_{l}\geq 0 of \operatorname{cone}(T,B) for ij\geq_{B}kl. For 4. notice that y satisfies y_{ij}=0 for each contracted edge ij if and only if x satisfies x_{i}=x_{j} when i and j label the same vertex in T(B); these are precisely the equalities that hold in \operatorname{cone}(T(B),B), as we explained in the proof of Lemma [2.9](https://arxiv.org/html/2603.03425v3#S2.Thmtheorem9 "Lemma 2.9. ‣ 2.3.3 The cosmohedral fan in terms of bracketed trees. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron"). ∎

###### Example 3.4.

Consider the Matryoshka M_{1} and its corresponding tree T_{1} with edges e=21,f=32,g=34, and bracketing B_{1}=\{e,efg,g\}, as shown in Figure [9](https://arxiv.org/html/2603.03425v3#S2.F9 "Figure 9 ‣ 2.3.3 The cosmohedral fan in terms of bracketed trees. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron"). Then

\operatorname{cone}_{T_{1}}(B_{1})_{\geq 0}=\{y\in\mathbb{R}^{efg}\,:\,y_{f}\geq\{y_{e},y_{g}\}\geq 0\}\subset\mathbb{R}^{E(T)},

and its image under the map f_{T_{1}} is

\operatorname{cone}(T_{1},B_{1})=\{x\in\mathbb{R}^{4}/\mathbb{R}\mathbf{1}\,:\,x_{3}-x_{2}\geq\{x_{2}-x_{1},x_{3}-x_{4}\}\geq 0\}\subset\mathbb{R}^{n}/\mathbb{R}\mathbf{1}

as illustrated in Figure [12](https://arxiv.org/html/2603.03425v3#S3.F12 "Figure 12 ‣ 3.2 The cosmohedral fan and positive bracket associahedral fans ‣ 3 Proof of correctness of the cosmohedral fan ‣ Combinatorics of the cosmohedron").

![Image 3: Refer to caption](https://arxiv.org/html/2603.03425v3/SubFanLabels.png)

Figure 12: The linear map f_{T_{1}} embeds the four-sided \operatorname{cone}_{T_{1}}(B_{1})_{\geq 0} of \operatorname{BrAssocFan}_{\geq 0}(T_{1}) in Figure [11](https://arxiv.org/html/2603.03425v3#S3.F11 "Figure 11 ‣ 3.1 Positive bracket associahedral fans ‣ 3 Proof of correctness of the cosmohedral fan ‣ Combinatorics of the cosmohedron") into the \operatorname{cone}(T_{1}) of \operatorname{Assoc}_{n-1} in Figure [6](https://arxiv.org/html/2603.03425v3#S2.F6 "Figure 6 ‣ 2.3.1 Preliminaries: plane trees, triangulations, and the associahedral fan. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron"). The result is the cone f_{T_{1}}(\operatorname{cone}_{T_{1}}(B_{1})_{\geq 0})=\operatorname{cone}(T_{1},B_{1}) of \operatorname{Cosmo}_{n-1}.

For the moment, we only know \operatorname{CosmoFan}_{n-1} to be the set of Matryoshkal cones. Our next goal is to show that it is indeed a polyhedral fan, with the correct combinatorial structure. We begin with a simple technical lemma.

###### Corollary 3.5.

We have

\operatorname{CosmoFan}_{n-1}=\bigcup_{T}f_{T}(\operatorname{BrAssocFan}_{\geq 0}(T))

where the union is over all binary trees T with n+1 leaves.

###### Proof.

Any cone of \operatorname{CosmoFan}_{n-1} is of the form \operatorname{cone}(M)=\operatorname{cone}(T,B) for a bracketed tree (T,B). If we choose any binary T^{+}\geq T, then B is a pre-bracketing of T^{+}, and Lemma [3.3](https://arxiv.org/html/2603.03425v3#S3.Thmtheorem3 "Lemma 3.3. ‣ 3.2 The cosmohedral fan and positive bracket associahedral fans ‣ 3 Proof of correctness of the cosmohedral fan ‣ Combinatorics of the cosmohedron").4 shows that this cone is in f_{T^{+}}(\operatorname{BrAssocFan}_{\geq 0}(T^{+})). ∎

Notice in particular that the union above is not disjoint.

### 3.3 The cosmohedral fan is a complete fan

###### Theorem 3.6.

The cosmohedral fan \operatorname{CosmoFan}_{n-1} is a complete fan in \mathbb{R}^{n}/\mathbb{R}\mathbf{1}. The map M\mapsto\operatorname{cone}(M) is an order-preserving bijection between the Matryoshkas of an (n+2)-gon and the faces of \operatorname{CosmoFan}_{n-1}.

###### Proof.

To prove that the cones intersect properly, let us first analyze how a Matryoshkal \operatorname{cone}(T_{1},B_{1}), which sits in an associahedral \operatorname{cone}(T_{1}) of the same dimension, intersects a smaller associahedral \operatorname{cone}(T)\subseteq\operatorname{cone}(T_{1}). Consider any maximal T^{+}\geq T_{1}. Then B_{1} is a pre-bracketing of T^{+} and we have

\displaystyle\operatorname{cone}(T_{1},B_{1})\cap\operatorname{cone}(T)\displaystyle=\displaystyle f_{T^{+}}(\operatorname{cone}_{T^{+}}(B_{1})_{\geq 0})\cap f_{T^{+}}(\mathbb{R}_{\geq 0}^{E(T)})
\displaystyle=\displaystyle f_{T^{+}}(\operatorname{cone}_{T^{+}}(B_{1})_{\geq 0}\cap\mathbb{R}_{\geq 0}^{E(T)})
\displaystyle=\displaystyle f_{T^{+}}(\operatorname{cone}_{T^{+}}(B_{1}/\overline{F})_{\geq 0})\quad\text{ for }F=E(T^{+})-E(T),\text{ by }\eqref{eq:brackettoprebracket}
\displaystyle=\displaystyle\operatorname{cone}(T^{\prime}_{1},B^{\prime}_{1})\quad\text{ for }B_{1}^{\prime}=B_{1}/\overline{F}

where T^{\prime}_{1}=T^{+}/\overline{F}=T/(\overline{F}-F)\leq T is the tree supporting B_{1}^{\prime}.

For example, in the top row of Figure [13](https://arxiv.org/html/2603.03425v3#S3.F13 "Figure 13 ‣ Proof. ‣ 3.3 The cosmohedral fan is a complete fan ‣ 3 Proof of correctness of the cosmohedral fan ‣ Combinatorics of the cosmohedron") we compute \operatorname{cone}(M_{1})\cap\operatorname{cone}(T) for the Matryoshka M_{1}=(T_{1},B_{1}) of Figure [9](https://arxiv.org/html/2603.03425v3#S2.F9 "Figure 9 ‣ 2.3.3 The cosmohedral fan in terms of bracketed trees. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron") and the tree T of the shown subdivision. Here T^{+}=T_{1}, F=E(T^{+})-E(T)=\{21\}, and the downset in \tau_{1} is \overline{F}=F=\{21\}. Thus we obtain (T_{1}^{\prime},B_{1}^{\prime}) by contracting edge 21 from (T_{1},B_{1}). The resulting inequalities are x_{3}-x_{21}\geq x_{3}-x_{4}\geq 0, where x_{21} denotes the common value of x_{2}=x_{1}. In the second row we compute \operatorname{cone}(M_{2})\cap\operatorname{cone}(T) similarly.

Figure 13: The intersection of the two Matryoshkal cones M_{1} and M_{2} described in Figure [9](https://arxiv.org/html/2603.03425v3#S2.F9 "Figure 9 ‣ 2.3.3 The cosmohedral fan in terms of bracketed trees. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron"). We begin by computing T=T_{1}\wedge T_{2} and use it to find T_{1}^{\prime},T_{2}^{\prime} and the inherited bracketings B_{1}^{\prime} and B_{2}^{\prime}. From there we find B^{\prime}=B_{1}^{\prime}\wedge B_{2}^{\prime}. Then M^{\prime} is the Matryoshka of T^{\prime} and B^{\prime}.

Now consider two Matryoshkas M_{1}=(T_{1},B_{1}) and M_{2}=(T_{2},B_{2}). Their cones intersect inside \operatorname{cone}(T_{1})\cap\operatorname{cone}(T_{2})=\operatorname{cone}(T) where T=T_{1}\wedge T_{2} in \operatorname{Tree}_{n-1}. Then

\displaystyle\operatorname{cone}(M_{1})\cap\operatorname{cone}(M_{2})\displaystyle=\displaystyle(\operatorname{cone}(T_{1},B_{1})\cap\operatorname{cone}(T))\cap(\operatorname{cone}(T_{2},B_{2})\cap\operatorname{cone}(T))
\displaystyle=\displaystyle\operatorname{cone}(T^{\prime}_{1},B^{\prime}_{1})\cap\operatorname{cone}(T^{\prime}_{2},B^{\prime}_{2})\text{ where }T_{1}^{\prime},T_{2}^{\prime}\leq T
\displaystyle=\displaystyle f_{T^{+}}(\operatorname{cone}_{T^{+}}(B^{\prime}_{1})_{\geq 0}\cap\operatorname{cone}_{T^{+}}(B^{\prime}_{2})_{\geq 0})\,\,\,\text{ for any maximal }T^{+}\geq T
\displaystyle=\displaystyle f_{T^{+}}(\operatorname{cone}_{T^{+}}(B^{\prime})_{\geq 0})\quad\text{ for $B^{\prime}=B^{\prime}_{1}\wedge B^{\prime}_{2}$, by Proposition \ref{prop:BrAssocFan}}
\displaystyle=\displaystyle\operatorname{cone}(M^{\prime})\text{ for the Matryoshka $M^{\prime}$ of }(T^{\prime},B^{\prime})

where T^{\prime}\leq T^{+} is the tree supporting B^{\prime}=B^{\prime}_{1}\wedge B^{\prime}_{2}.

In the third row of Figure [13](https://arxiv.org/html/2603.03425v3#S3.F13 "Figure 13 ‣ Proof. ‣ 3.3 The cosmohedral fan is a complete fan ‣ 3 Proof of correctness of the cosmohedral fan ‣ Combinatorics of the cosmohedron") we compute \operatorname{cone}(M_{1})\cap\operatorname{cone}(M_{2}) for the Matryoshkas of Figure [9](https://arxiv.org/html/2603.03425v3#S2.F9 "Figure 9 ‣ 2.3.3 The cosmohedral fan in terms of bracketed trees. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron"), noticing that the tree we had considered is T=T_{1}\wedge T_{2}. Considering B_{1}^{\prime} and B_{2}^{\prime} as pre-bracketings of any maximal tree T^{+}\geq T we compute their meet B^{\prime}=B_{1}^{\prime}\wedge B_{2}^{\prime} and obtain from it the tree T^{\prime} and the inequalities of the intersection, x_{3}-x_{12}=x_{3}-x_{4}\geq 0; that is, the ray x_{3}\geq x_{124} in \mathbb{R}^{4}/\mathbb{R}.

To prove that the fan is complete, we begin with the fact that the maximal tree cones \operatorname{cone}(T) of the associahedral fan cover \mathbb{R}^{n}/\mathbb{R}\mathbf{1}. Now, for each T, the cones \operatorname{cone}_{T}(B)_{\geq 0} of the pre-bracketing fan cover the positive orthant \mathbb{R}^{E(T)}_{\geq 0}, so their images f_{T}(\operatorname{cone}_{T}(B)_{\geq 0})=\operatorname{cone}(T(B),B) cover f_{T}(\mathbb{R}^{E(T)}_{\geq 0})=\operatorname{cone}(T).

Finally, let us prove the poset isomorphism, keeping in mind that we already showed that \operatorname{Mat}_{n+2}\cong\operatorname{BrTree}_{n+1}. Clearly M_{1}\subseteq M_{2} implies \operatorname{cone}(M_{1})\subseteq\operatorname{cone}(M_{2}). Assume that \operatorname{cone}(M_{1})=\operatorname{cone}(T_{1},B_{1})\subseteq\operatorname{cone}(T_{2},B_{2})=\operatorname{cone}(M_{2}). Let T be a maximal tree such that \operatorname{cone}(T) contains these cones; B_{1} and B_{2} are pre-bracketings of T. Then f_{T}(\operatorname{cone}_{T}(B_{1})_{\geq 0})=\operatorname{cone}(T_{1},B_{1})=\operatorname{cone}(T_{2},B_{2})\subseteq f_{T}(\operatorname{cone}_{T}(B_{2})_{\geq 0}) by Lemma [3.3](https://arxiv.org/html/2603.03425v3#S3.Thmtheorem3 "Lemma 3.3. ‣ 3.2 The cosmohedral fan and positive bracket associahedral fans ‣ 3 Proof of correctness of the cosmohedral fan ‣ Combinatorics of the cosmohedron") so \operatorname{cone}_{T}(B_{1})_{\geq 0}\subseteq\operatorname{cone}_{T}(B_{2})_{\geq 0} which means that B_{1}\leq B_{2} in \operatorname{PreBr}(T) by Proposition [3.2](https://arxiv.org/html/2603.03425v3#S3.Thmtheorem2 "Proposition 3.2. ‣ 3.1 Positive bracket associahedral fans ‣ 3 Proof of correctness of the cosmohedral fan ‣ Combinatorics of the cosmohedron"). Therefore (T_{1},B_{1})\subseteq(T_{2},B_{2}) in \operatorname{BrTree}_{n-1}, which means that M_{1}\leq M_{2} in \operatorname{Mat}_{n-1} by Proposition [2.13](https://arxiv.org/html/2603.03425v3#S2.Thmtheorem13 "Proposition 2.13. ‣ 2.3.3 The cosmohedral fan in terms of bracketed trees. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron"). ∎

## 4 The cosmohedron

In this section, we define the cosmohedron \operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}) and prove that its normal fan is the cosmohedral fan \operatorname{CosmoFan}_{n-1}.

To set conventions, let e_{1},\ldots,e_{n} be the standard basis of \mathbb{R}^{n}. Let us write i<<j, when integers i<j are not consecutive, that is, |j-i|>1. We note that the set of pairs (i,j) with 0\leq i<<j\leq n+1 is in bijection with the set \operatorname{diag}^{+}(\pentagon_{n+2})=\operatorname{diag}(\pentagon_{n+2})\cup\{e\} of diagonals augmented with the base edge e.

### 4.1 Preliminaries: Associahedra and bracket associahedra

Loday’s associahedron. Given a sequence of positive numbers \mathsf{a}=(a_{ij})_{0\leq i<<j\leq n+1}, the _Loday associahedron_ is

\displaystyle\operatorname{Loday}_{n-1}(\mathsf{a})\displaystyle=\sum_{0\leq i<<j\leq n+1}a_{ij}\Delta_{(i,j)}(3)
\displaystyle=\{x\in\mathbb{R}^{n}:\sum_{1\leq i\leq n}x_{i}=\sum_{0\leq i<<j\leq n+1}a_{ij},
\displaystyle\qquad\qquad\quad\sum_{r\leq i\leq s}x_{i}\geq\sum_{r-1\leq i<<j\leq s+1}a_{ij}\text{ for }1\leq r\leq s\leq n.\}(4)

where \Delta_{(i,j)}=\Delta_{[i+1,j-1]}=\operatorname{conv}(e_{i+1},\ldots,e_{j-1}) is a face of the standard simplex in \mathbb{R}^{n}, and P+Q=\{p+q\,:\,p\in P,q\in Q\} denotes Minkowski sum of polytopes.1 1 1 It is more customary to choose \mathbf{\alpha}=(\alpha_{ij})_{1\leq i\leq j\leq n} and write \sum_{i\leq j}\alpha_{ij}\Delta_{[i,j]} for Loday’s associahedron; we use a_{i-1,j+1}=\alpha_{ij} instead.

###### Theorem 4.1.

[[AA23](https://arxiv.org/html/2603.03425v3#bib.bibx1), [Pos09](https://arxiv.org/html/2603.03425v3#bib.bibx32)] The normal fan of the associahedron \operatorname{Loday}_{n-1}(\mathsf{a}) is the associahedral fan \operatorname{AssocFan}_{n-1}.

Let \mathcal{S} be a triangulation of \pentagon_{n+2} corresponding to the binary tree T. It will be useful to express the vertex a_{T} of Loday’s associahedron in terms of \mathcal{S}. To do so, say that diagonal ij _transversally crosses_ the triangle labeled v with vertices u<v<w if it intersects edges uv and vw, not at v. This is equivalent to u\leq i<v<j\leq w. When this is the case we write ij\rightarrow v. We invite the reader to check that each diagonal of the polygon crosses exactly one triangle of the triangulation \mathcal{S} transversally.

###### Lemma 4.2.

[[Lod04](https://arxiv.org/html/2603.03425v3#bib.bibx29), [Pos09](https://arxiv.org/html/2603.03425v3#bib.bibx32)] The vertex a_{\mathcal{S}} of \operatorname{Loday}_{n-1}(\mathsf{a}) corresponding to a triangulation \mathcal{S} of \pentagon_{n+2} has coordinates

(a_{\mathcal{S}})_{v}=\sum_{ij\rightarrow v\text{ in }\mathcal{S}}a_{ij}

summing over the extended diagonals ij\in\operatorname{diag}^{+}(\pentagon_{n+2}) that cross triangle v transversally.

###### Proof.

The vertex a_{T} of \operatorname{Loday}_{n-1}(\mathsf{a}) corresponding to a binary tree T has coordinates

(a_{T})_{v}=\sum_{i\vee j=v\text{ in }T}a_{i-1,j+1}\qquad\text{ for }v\in[n]

where i\vee j is the lowest common ancestor of i and j; this is implicit in [[Lod04](https://arxiv.org/html/2603.03425v3#bib.bibx29), [Pos09](https://arxiv.org/html/2603.03425v3#bib.bibx32)] see for example [[AA23](https://arxiv.org/html/2603.03425v3#bib.bibx1), Figure 6]. The second statement then follows readily from the bijection. ∎

Devadoss’s bracket associahedron. For a graph H, let \mathsf{b}:\operatorname{Brackets}(H)\rightarrow\mathbb{Z} be any function that is _concave_ in the sense that b(\beta_{1})+b(\beta_{2})>b(\beta_{1}\cap\beta_{2})+b(\beta_{1}\cup\beta_{2}) for any intersecting brackets \beta_{1} and \beta_{2}. For example, b(\beta) can be any concave function of the number of edges |\beta|=k; that is, b(k)>\frac{1}{2}(b(k-1)+b(k+1)).

###### Definition 4.3.

The _Devadoss bracket associahedron_ of a graph H is

\displaystyle\operatorname{BrAssoc}_{H}(\mathsf{b})=\{y\in\mathbb{R}^{E(H)}\displaystyle:\sum_{e\in E(H)}y_{e}=b(H),
\displaystyle\,\,\,\,\,\,\,\sum_{e\in\beta}y_{e}\leq b(\beta)\text{ for all brackets $\beta$ of }H\}.(5)

###### Theorem 4.4.

[[Dev09](https://arxiv.org/html/2603.03425v3#bib.bibx20), [PPP23](https://arxiv.org/html/2603.03425v3#bib.bibx33)] The normal fan of the bracket associahedron \operatorname{BrAssoc}_{H}{\mathsf{b}} is the bracket associahedral fan of H.

We remark that Devadoss [[Dev09](https://arxiv.org/html/2603.03425v3#bib.bibx20)] works with the sequence b(\beta)=-3^{|\beta|-2} – see Lemma [5.2](https://arxiv.org/html/2603.03425v3#S5.Thmtheorem2 "Lemma 5.2. ‣ 5 Face structure and enumeration ‣ Combinatorics of the cosmohedron") – but in fact Theorem [4.4](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem4 "Theorem 4.4. ‣ 4.1 Preliminaries: Associahedra and bracket associahedra ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron") holds precisely for concave sequences. To see this, we use Padrol, Pilaud, and Poullot’s description [[PPP23](https://arxiv.org/html/2603.03425v3#bib.bibx33), Corollary 2.12] of the deformation cone of graph associahedral fans, and in particular bracket associahedral fans. This cone is given by the equality b(G)=0 and inequalities b(\beta-e)+b(\beta-e^{\prime})<b(\beta)+b(\beta-e-e^{\prime}) for any bracket \beta and any non-disconnecting edges e,e^{\prime} of \beta. This is equivalent to the concavity condition above.2 2 2 Padrol, Pilaud, and Poullot require b(G)=0, but we can remove this condition by a suitable translation. For general graphs, some care is required when \beta-e-e^{\prime} is disconnected, but this does not happen for trees.

###### Lemma 4.5.

Let \mathsf{b}:\operatorname{Brackets}(H)\rightarrow\mathbb{R} be concave. The vertex b_{B} of \operatorname{BrAssoc}_{H}(\mathsf{b}) associated to a maximal bracketing B of H is given by

(b_{B})_{e}=b(\operatorname{bracket}(e))-\sum_{f\lessdot e}b(\operatorname{bracket}(f))\qquad\text{ for }e\in E(H),(6)

summing over the edges f that e covers in \operatorname{tree}(B).

In a maximal bracketing, each bracket contains one or two maximal sub-brackets, so the expression above has at most three terms.

###### Proof.

For each maximal bracketing B, we can compute the entries of vertex b_{B} recursively by (b_{B})_{e}=b(e) for any edge e that is a bracket, and \sum_{e\in\beta}(b_{B})_{e}=b(\beta) for every bracket \beta of B. In the tree of B, this says b(\operatorname{bracket}(e))=\sum_{f\leq e}(b_{B})_{e}. Adding a minimum element \widehat{0} with b(\widehat{0})=b_{B}(\widehat{0})=0, and applying Möbius inversion in the poset \operatorname{tree}(B)\cup\widehat{0}, we get (b_{B})_{e}=\sum_{f\leq e}\mu(f,e)\,b(\operatorname{bracket}(f)) Since the interval [f,e] is a path, \mu(f,e) equals 1 if f=e, it equals -1 if f\lessdot e, and it is 0 otherwise; so ([6](https://arxiv.org/html/2603.03425v3#S4.E6 "In Lemma 4.5. ‣ 4.1 Preliminaries: Associahedra and bracket associahedra ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron")) does indeed describe Devadoss’s bracket associahedron. ∎

Slightly differently, we will work with the _extended set of brackets_\text{Brackets}^{+}(H)=\text{Brackets}(H)\cup V(H), including also a _small bracket_ around each individual vertex. Consider any concave function \mathsf{b}:\text{Brackets}^{+}(T)\rightarrow\mathbb{Z}; in the concavity condition, if \beta_{1}\cap\beta_{2} is an edge-less graph on vertices v_{1},\ldots,v_{k}, we set b(\beta_{1}\cap\beta_{2})=b(v_{1})+\cdots+\beta(v_{k}). The function \mathsf{b} induces a concave function \mathsf{b}^{-}:\operatorname{Brackets}(H)\rightarrow\mathbb{Z} defined by

b^{-}(\beta)=b(\beta)-\sum_{v\in V(\beta)}b(v),(7)

whose bracket associahedron we now describe. A bracketing B of H now gives an _extended_\operatorname{tree}^{+}(B) that also includes the small brackets.

###### Lemma 4.6.

Let \mathsf{b}:\operatorname{Brackets}^{+}(H)\rightarrow\mathbb{R} be concave and \mathsf{b}^{-}:\operatorname{Brackets}(H)\rightarrow\mathbb{R} be given by ([7](https://arxiv.org/html/2603.03425v3#S4.E7 "In 4.1 Preliminaries: Associahedra and bracket associahedra ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron")). The vertex b_{B} of \operatorname{BrAssoc}_{H}(\mathsf{b}^{-}) associated to a maximal bracketing B of H is given by

(b_{B})_{e}=b(\operatorname{bracket}(e))-\sum_{f\lessdot e}b(\operatorname{bracket}(f))\qquad\text{ for }e\in E(H).(8)

summing over the edges f that e covers in the extended \operatorname{tree}^{+}(B).

###### Proof.

This follows by applying Lemma [4.5](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem5 "Lemma 4.5. ‣ 4.1 Preliminaries: Associahedra and bracket associahedra ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron") to \mathsf{b}^{-}. ∎

### 4.2 The cosmohedron

In this section we construct a _cosmohedron_: a polytope \operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}) whose normal fan is the cosmohedral fan \operatorname{CosmoFan}_{n-1}. We first give two vertex descriptions of \operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}) and use them to prove that it has the correct normal fan. We then use the normal fan to compute the inequality description of \operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}). Finally, we use those inequalities to prove that \operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}) is linearly isomorphic to the cosmohedron proposed by Arkani-Hamed, Figueiredo, and Vazão [[AHFV25](https://arxiv.org/html/2603.03425v3#bib.bibx8)]. This gives a proof of the correctness of their construction.

Recall that \operatorname{diag}^{+}(\pentagon_{n+2})=\operatorname{diag}(\pentagon_{n+2})\cup e is the set of diagonals of the (n+2)-gon augmented with the base e. We say a function \mathsf{b}=(b(P)\,:\,P\text{ is a subpolygon of }\pentagon_{{n+2}})\in\mathbb{R}^{\binom{n+2}{\geq 3}} is _concave_ if b(P_{1})+b(P_{2})>b(P_{1}\cap P_{2})+b(P_{1}\cup P_{2}) whenever the intersection and union are also subpolygons. An example is a concave function of the number of vertices of the polygon.

Throughout this section, we fix:

\bullet an arbitrary positive vector \mathsf{a}=(a_{ij}\,:\,ij\in\operatorname{diag}^{+}(\pentagon_{n+2}))\in\mathbb{R}^{\operatorname{diag}^{+}(\pentagon_{n+2})}\cong\mathbb{R}^{\binom{n+1}{2}},

\bullet a concave vector \mathsf{b}=(b(P)\,:\,P\text{ is a subpolygon of }\pentagon_{{n+2}})\in\mathbb{R}^{\binom{n+2}{\geq 3}},

\bullet a small constant 0<\varepsilon<\frac{m}{24M} where m=\min a_{ij} and M=\max|b(P)|.

The cosmohedron \operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}) is constructed from those choices.

#### 4.2.1 Vertex description

For a rooted plane tree T let

\displaystyle f^{T}\,:\,\mathbb{R}^{E(T)}\displaystyle\longrightarrow\displaystyle\mathbb{R}^{V(T)}_{0}
\displaystyle e_{ij}\displaystyle\longmapsto\displaystyle e_{i}-e_{j}

be the isomorphism between \mathbb{R}^{E(T)} and \mathbb{R}^{V(T)}_{0}:=\{x\in\mathbb{R}^{V(T)}\,:\,\sum_{v}x_{v}=0\} dual to the isomorphism f_{T} of Lemma [3.3](https://arxiv.org/html/2603.03425v3#S3.Thmtheorem3 "Lemma 3.3. ‣ 3.2 The cosmohedral fan and positive bracket associahedral fans ‣ 3 Proof of correctness of the cosmohedral fan ‣ Combinatorics of the cosmohedron"). Also, our chosen vector \mathsf{b}\in\mathbb{R}^{\binom{n+2}{\geq 3}} induces a map \mathsf{b}:\operatorname{Brackets}^{+}(T)\rightarrow\mathbb{R}, since the extended brackets of T correspond to the polygons in the triangulation of M. In turn this induces a map \mathsf{b}^{-}:\operatorname{Brackets}(T)\rightarrow\mathbb{R} by ([7](https://arxiv.org/html/2603.03425v3#S4.E7 "In 4.1 Preliminaries: Associahedra and bracket associahedra ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron")).

###### Definition 4.7.

For a maximal Matryoshka M with bracketed tree (T,B) we define

c_{M}:=a_{T}+\varepsilon\,f^{T}(b_{B})

where a_{T} is the vertex of \operatorname{Loday}_{n-1}(\mathsf{a}) corresponding to the plane tree T and b_{B} is the vertex of \operatorname{BrAssoc}_{T}(\mathsf{b}^{-}) corresponding to its bracketing B. The _cosmohedron_ is

\operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}):=\operatorname{conv}(c_{M}\,:\,M\text{ is a maximal Matryoshka on the $(n+2)$-gon})\subset\mathbb{R}^{n}.

Before we prove the correctness of this construction, let us give an explicit formula for the vertex c_{M} directly in terms of the Matryoshka M.

###### Proposition 4.8.

Let M be a maximal Matryoshka on \pentagon_{n+2}. For each v\in[n] define

(a_{M})_{v}=\sum_{ij\rightarrow v}a_{ij}=\sum_{u\leq i<v}\sum_{v<j\leq w}a_{ij}

summing over the diagonals ij that cross triangle v (with vertices uvw) of the triangulation of M transversally.

For each polygon of the Matryoshka M that is subdivided into maximal subpolygons Q and R in M, let

b_{M}(P)=b(P)-b(Q)-b(R).

For each v\in[n], where triangle v has edges u_{v} (opposite to vertex v), l_{v}, and r_{v}, define

(b_{M})_{v}=b_{M}(P(l_{v}))+b_{M}(P(r_{v}))-b_{M}(P(u_{v}))

where P(e) is the minimal polygon in M containing diagonal e; when e is an edge of \pentagon_{n+2}, we set b_{M}(P(e))=0. Then the vertex c_{M} of a Matryoshka M has coordinates

(c_{M})_{v}=(a_{M})_{v}+\varepsilon(b_{M})_{v}\qquad\text{ for }v\in[n].

###### Proof.

This follows readily from Lemmas [4.2](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem2 "Lemma 4.2. ‣ 4.1 Preliminaries: Associahedra and bracket associahedra ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron") and [4.5](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem5 "Lemma 4.5. ‣ 4.1 Preliminaries: Associahedra and bracket associahedra ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron") and Definition [4.7](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem7 "Definition 4.7. ‣ 4.2.1 Vertex description ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron"). ∎

###### Example 4.10.

Let us compute the vertex c_{M_{2}} for the Matryoshka M_{2} of Figure [4](https://arxiv.org/html/2603.03425v3#S2.F4 "Figure 4 ‣ 2.1 Matryoshkas ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron"). We have a_{M_{2}}=(1\cdot 2,1\cdot 1,3\cdot 2,1\cdot 1)=(2,1,6,1) The non-zero coordinates of b_{M_{2}} are b_{M_{2}}(\text{hexagon})=3-1-0=2 and b_{M_{2}}(\text{pentagon})=1-0-0=1. Therefore b_{M_{2}}=(0+0-1,0+0-0,1+2-0,0+0-2)=(-1,0,3,-2) and

c_{M_{2}}=(2,1,6,1)+\varepsilon(-1,0,3,2).

#### 4.2.2 The cosmohedron: proof of correctness

###### Theorem 4.11.

The normal fan of the cosmohedron \operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}) is the cosmohedral fan \operatorname{CosmoFan}_{n-1}.

###### Proof.

Let \mathcal{N}_{P}(v):=\{w\in V^{*}\,:\,w(v)\geq w(p)\text{ for all }p\in P\} denote the normal cone of vertex v in polytope P\subset V. Since the cosmohedron lives on a hyperplane with constant coordinate sum, its dual space is \mathbb{R}^{n}/\mathbb{R}\mathbf{1}. We need to prove that

\mathcal{N}_{\operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b})}(c_{M})=\operatorname{cone}(M)

for all maximal Matryoshkas M. Since both the left-hand sides and the right-hand sides subdivide \mathbb{R}^{n}/\mathbb{R}\mathbf{1} with no interior overlap as M ranges over the maximal Matryoshkas, it suffices to prove the inclusion \subseteq for all M.

Let M be a maximal Matryoshka with bracketed tree (T,B), and take w\in\mathcal{N}_{\operatorname{Cosmo}_{n-1}}(c_{M}). To prove that w\in\operatorname{cone}(M), we proceed in two steps.

Step 1. Consider a covering relation \pentagon_{M}(ij^{\perp})\gtrdot\pentagon_{M}(kl^{\perp}) in M. We need to show that w satisfies the corresponding inequality w_{i}-w_{j}\geq w_{k}-w_{l} for \operatorname{cone}(M).

Edge kl^{\perp} creates a subdivision \pentagon_{M}(kl^{\perp})=P_{1}\cup P_{2} and edge ij^{\perp} creates a subdivision \pentagon_{M}(ij^{\perp})=\pentagon_{M}(kl^{\perp})\cup P_{3}. Consider the Matryoshka M^{\prime} obtained by removing \pentagon_{M}(kl^{\perp}) and adding a new polygon around P_{2} and P_{3}.

Figure 14: A swap of a non-minimal polygon in a Matryoshka, and the corresponding effect on the polygon containment tree \tau. 

The bracketed tree of M^{\prime} is (T,B^{\prime}) for the same tree T and a different bracketing. The containment trees \tau=\operatorname{tree}(B) and \tau^{\prime}=\operatorname{tree}(B^{\prime}) differ by a simple mutation, as shown in Figure [14](https://arxiv.org/html/2603.03425v3#S4.F14 "Figure 14 ‣ Proof. ‣ 4.2.2 The cosmohedron: proof of correctness ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron"). Since Devadoss’s bracket associahedron is a generalized permutahedron, its adjacent vertices b_{B} and b_{B^{\prime}} satisfy b_{B}-b_{B}^{\prime}=L(e_{ij}-e_{kl}) for some L>0[[AA23](https://arxiv.org/html/2603.03425v3#bib.bibx1), [Pos09](https://arxiv.org/html/2603.03425v3#bib.bibx32)]. Now w\in\mathcal{N}_{\operatorname{Cosmo}_{n-1}}(c_{M}) gives w(c_{M})\geq w(c_{M^{\prime}}), so

\displaystyle w(a_{T}+\varepsilon\,f^{T}(b_{B}))\displaystyle\geq\displaystyle w(a_{T}+\varepsilon\,f^{T}(b_{B^{\prime}}))
\displaystyle wf^{T}(L(e_{ij}-e_{kl}))\displaystyle\geq\displaystyle 0
\displaystyle w((e_{i}-e_{j})-(e_{k}-e_{l}))\displaystyle\geq\displaystyle 0

as desired.

Step 2. Consider a polygon \pentagon_{M}(kl^{\perp}) that does not contain a smaller \pentagon_{M}(k^{\prime}l^{\prime\perp}); it must be a quadrilateral split into two triangles k and l. Consider the Matryoshka M^{\prime} obtained by flipping those two triangles in the quadrilateral.

Figure 15: A flip of the triangulation in a Matryoshka, and the corresponding effect on the tree T of the triangulation. 

The tree T^{\prime} of (the triangulation of) M^{\prime} is obtained from T by a simple edge mutation, as shown in Figure [15](https://arxiv.org/html/2603.03425v3#S4.F15 "Figure 15 ‣ Proof. ‣ 4.2.2 The cosmohedron: proof of correctness ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron"). Since Loday’s associahedron is a generalized permutahedron, its adjacent vertices a_{T} and a_{T^{\prime}} satisfy a_{T}-a_{T}^{\prime}=L(e_{i}-e_{j}) for some L>0[[AA23](https://arxiv.org/html/2603.03425v3#bib.bibx1), [Pos09](https://arxiv.org/html/2603.03425v3#bib.bibx32)]. The bracketings B and B^{\prime} both contain bracket ij, so they contain the same set of polygons, and give isomorphic posets \tau=\operatorname{tree}(B)\cong\operatorname{tree}(B^{\prime})=\tau^{\prime}. This means that b_{B}=b_{B^{\prime}} up to the resulting reordering of the edges (hi,ij,jk)\mapsto(hj,ji,ik). Therefore if we write H=(b_{B})_{hi},I=(b_{B})_{ij},J=(b_{B})_{jk},

\displaystyle f^{T}(b_{B})-f^{T^{\prime}}(b_{B^{\prime}})\displaystyle=\displaystyle(H(e_{h}-e_{i})+I(e_{i}-e_{j})+J(e_{j}-e_{k}))
\displaystyle\quad-(H(e_{h}-e_{j})+I(e_{j}-e_{i})+J(e_{i}-e_{k}))
\displaystyle=\displaystyle(-H+2I-J)(e_{i}-e_{j}).

Now since w\in\mathcal{N}_{\operatorname{Cosmo}_{n-1}}(c_{M}), we have w(c_{M})\geq w(c_{M^{\prime}}), so

\displaystyle w\left((a_{T}-a_{T^{\prime}})+\varepsilon\,(f^{T}-f^{T^{\prime}})(b_{B})\right)\displaystyle\geq\displaystyle 0
\displaystyle w((L+2\varepsilon(-H+2I-J))(e_{i}-e_{j})\displaystyle\geq\displaystyle 0
\displaystyle w(e_{i}-e_{j})\displaystyle\geq\displaystyle 0

as desired. In the last step we are observing that each entry of (b_{B}) is at most 3M because the formula in Lemma [4.5](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem5 "Lemma 4.5. ‣ 4.1 Preliminaries: Associahedra and bracket associahedra ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron") contains at most three terms, so |2\varepsilon(-H+2I-J)|<24M\varepsilon, whereas each edge a_{T}-a_{T^{\prime}} of \operatorname{Loday}_{n-1}(\mathsf{a}) is a sum of edges of its Minkowski summands, so its length L is at least m. ∎

#### 4.2.3 Inequality description.

For a chord C and a vertex i or diagonal ij of \pentagon_{n+2}, let us write C\succ i if the chord C is strictly above vertex i and C\succeq ij if the chord C is weakly above diagonal ij, with respect to the base edge e at the top.

For a subdivision \mathcal{S} of \pentagon_{n+2} (or equivalently a non-crossing set of chords \mathcal{C}), we define some parameters of \mathcal{S} and \mathcal{C} in a few equivalent ways. Each description is sufficient, but it can be useful at times to have all of them available. Let the _depth_ of a vertex i or a diagonal ij of \pentagon_{n+2} be

\displaystyle d_{\mathcal{S}}(i)\displaystyle=\displaystyle\#\text{ of chords $C$ of $\mathcal{C}$ with $C\succ i$}
\displaystyle=distance from vertex i to the root in \operatorname{tree}(\mathcal{S})
\displaystyle=\displaystyle\text{depth of vertex $i$ in $\operatorname{tree}_{\square}(\mathcal{S})$},\text{ and}
\displaystyle d_{\mathcal{S}}(ij)\displaystyle=\displaystyle\#\text{ of chords $C$ of $\mathcal{C}$ with $C\succeq ij$}
\displaystyle=\displaystyle\text{depth of cell $ij$ in $\operatorname{tree}_{\square}(\mathcal{S})$}.

where \operatorname{tree}_{\square}(\mathcal{S}) is the canonical embedding of the tree of S in a square grid, and the depth of a vertex or a cell in \operatorname{tree}_{\square}(\mathcal{S}) is the number of non-root vertices of the tree whose southwest/southeast shadow contains it; see Figure [16](https://arxiv.org/html/2603.03425v3#S4.F16 "Figure 16 ‣ 4.2.3 Inequality description. ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron").

Figure 16: A subdivision \mathcal{S}, its corresponding tree, its grid embedding with the corresponding cell depths, and the corresponding cosmohedral facet.

###### Theorem 4.13.

The cosmohedron \operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b})\subset\mathbb{R}^{n} is described by the irredundant system of inequalities

\displaystyle\sum_{1\leq i\leq n}x_{i}\displaystyle=\displaystyle\sum_{ij\in\operatorname{diag}(\pentagon_{n+2})}a_{ij},
\displaystyle\sum_{1\leq i\leq n}d_{\mathcal{S}}(i)\ x_{i}\displaystyle\geq\displaystyle\sum_{ij\in\operatorname{diag}(\pentagon_{n+2})}d_{\mathcal{S}}(ij)\ a_{ij}-\varepsilon\,\sum_{P\in\mathcal{S}}b(P)\text{ for any subdivision }\mathcal{S}.

###### Proof.

The lineality space of \operatorname{CosmoFan}_{n-1} is \mathbb{R}(1,\ldots,1), so the only equality satisfied by \operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}) is the given one.

The irredundant inequalities of \operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}) correspond to the rays of \operatorname{CosmoFan}_{n-1}, which are in bijection with the minimal non-trivial Matryoshkas by Theorem [3.6](https://arxiv.org/html/2603.03425v3#S3.Thmtheorem6 "Theorem 3.6. ‣ 3.3 The cosmohedral fan is a complete fan ‣ 3 Proof of correctness of the cosmohedral fan ‣ Combinatorics of the cosmohedron"). These are the subdivisions \mathcal{S} of \pentagon_{n+2} – with no nesting – obtained from the non-crossing sets of chords \mathcal{C}. By ([1](https://arxiv.org/html/2603.03425v3#S2.E1 "In 2.2 The cosmohedral fan in terms of Matryoshkas ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron")), the ray \mathbb{R}_{\geq 0}(\rho_{\mathcal{S}})=\operatorname{cone}(\mathcal{S}) is given by the inequalities x_{v}-x_{w}\geq 0, and hence (\rho_{\mathcal{S}})_{v}-(\rho_{\mathcal{S}})_{w}=1, for all chords vw^{\perp} of \mathcal{S}. Therefore \rho_{\mathcal{S}}=(d_{\mathcal{S}}(1),\ldots,d_{\mathcal{S}}(n)). Notice that \displaystyle\rho_{\mathcal{S}}=\sum_{C\in\mathcal{C}}{\rho_{\{C\}}} and \rho_{\{C\}}=e_{i+1}+\cdots+e_{j-1} for C=ij.

To determine the extreme value of (\rho_{\mathcal{S}},x) for x\in\operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}), we need to compute (\rho_{\mathcal{S}},c_{M}) for any maximal Matryoshka M=(T,B) that refines \mathcal{S}; see Figure [17](https://arxiv.org/html/2603.03425v3#S4.F17 "Figure 17 ‣ 4.2.3 Inequality description. ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron"). A chord C gets a label vw^{\perp} for a pair of adjacent triangles v,w in the triangulation T, and the definition in Lemma [3.3](https://arxiv.org/html/2603.03425v3#S3.Thmtheorem3 "Lemma 3.3. ‣ 3.2 The cosmohedral fan and positive bracket associahedral fans ‣ 3 Proof of correctness of the cosmohedral fan ‣ Combinatorics of the cosmohedron") gives f_{T}(e_{vw})=\rho_{\{C\}}+\mathbb{R}\mathbf{1}; this is also a ray of the associahedral fan \operatorname{AssocFan}_{n-1}. Therefore

\displaystyle(\rho_{\mathcal{C}},c_{M})\displaystyle=\displaystyle\sum_{C\in\mathcal{C}}\left(\rho_{C}\,,\,a_{T}+\varepsilon\,f^{T}(b_{B})\right)\quad\text{ by Definition \ref{def:vertmat},}
\displaystyle=\displaystyle\sum_{C\in\mathcal{C}}(\rho_{C},a_{T})+\varepsilon\,\sum_{C=vw^{\perp}\in\mathcal{C}}(f_{T}(e_{vw}),f^{T}(b_{B}))

The first sum comes from the Loday associahedron; using ([4](https://arxiv.org/html/2603.03425v3#S4.E4 "In 4.1 Preliminaries: Associahedra and bracket associahedra ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron")) we obtain

\sum_{C\in\mathcal{C}}(\rho_{C},a_{T})=\sum_{C}\sum_{ij\,:\,ij\preceq C}a_{ij}\\
=\sum_{ij}d_{\mathcal{C}}(ij)a_{ij}

Since f_{T} and f^{T} are dual, the second sum equals

\sum_{vw^{\perp}\in\mathcal{C}}(b_{B})_{vw}=\sum_{vw^{\perp}\in\mathcal{C}}\left(b(\operatorname{bracket}_{B}(vw))-\sum_{tu\lessdot_{B}vw}b(\operatorname{bracket}_{B}(tu))\right)

using Lemma [4.6](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem6 "Lemma 4.6. ‣ 4.1 Preliminaries: Associahedra and bracket associahedra ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron") for the vertices of the bracket associahedron. Here tu\lessdot_{B}vw means that vw is covers tu in the extended \operatorname{tree}^{+}(B) of the bracketing B.

Now, every diagonal d of the triangulation T that is not in \mathcal{C} is inside a polygon P\in\mathcal{S}, so \operatorname{bracket}_{B}(d)\subset P cannot contain any chord in \mathcal{C}. It follows that \mathcal{C} is an upper ideal of \operatorname{tree}^{+}(B), and that the brackets of the maximal non-elements of \mathcal{C} are precisely the polygons in \mathcal{S}. Also, \operatorname{tree}(B) has a unique maximal element VW=\operatorname{tree}(B)_{max}. Therefore

\displaystyle\sum_{vw^{\perp}\in\mathcal{C}}(b_{B})_{vw}\displaystyle=\displaystyle b(\operatorname{bracket}_{B}(VW))-\sum_{tu\in(\operatorname{tree}(B)-\mathcal{C})_{\text{max}}}b(\operatorname{bracket}_{B}(tu))
\displaystyle=\displaystyle b(\pentagon)-\sum_{P\in\mathcal{S}}b(P)
\displaystyle=\displaystyle-\sum_{P\in\mathcal{S}}b(P)

using that b(\pentagon)=0 . The desired result follows. ∎

Figure 17: A vertex on the facet given by Figure [16](https://arxiv.org/html/2603.03425v3#S4.F16 "Figure 16 ‣ 4.2.3 Inequality description. ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron"). 

###### Example 4.14.

Figure [17](https://arxiv.org/html/2603.03425v3#S4.F17 "Figure 17 ‣ 4.2.3 Inequality description. ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron") shows a Matryoshka M that refines the subdivision \mathcal{S} of Figure [16](https://arxiv.org/html/2603.03425v3#S4.F16 "Figure 16 ‣ 4.2.3 Inequality description. ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron"), so vertex c_{M} is on the facet that maximizes \rho_{\mathcal{C}}. In the computation of (\rho_{\mathcal{C}},c_{M}) above, the left sum is simply the dot product of the cell labels of \operatorname{tree}_{\square}(\mathcal{S}) and \mathsf{a}. Writing b(\cdot)=b(\operatorname{bracket}_{B}(\cdot)) for legibility, the right sum is

\displaystyle[b(42))-b(23)-b(45)]+[b(23)-b(3)-b(21)]+[b(45)-b(56)-b(4)]
\displaystyle=\displaystyle 0-b(3)-b(21)-b(56)-b(4)
\displaystyle=\displaystyle-b_{0124}-b_{234}-b_{047}-b_{4567}

#### 4.2.4 Isomorphism with the Arkani–Hamed-Figueiredo-Vazão cosmohedron.

In _kinematic space_ one assigns a real number to each diagonal of \pentagon_{n+2}:

\mathbb{R}^{\operatorname{diag}(\pentagon_{n+2})}=\{(X_{ij})_{ij\in\operatorname{diag}(\pentagon_{n+2})}\}\cong\mathbb{R}^{(n+2)(n-1)/2}

We assign X_{0,n+1}=0 to the base edge e=0,n+1.

Let us fix:

\bullet an arbitrary assignment \mathsf{A}=(A_{ij}\,:\,0\leq i<<j\leq n) of {\binom{n}{2}} positive numbers to the diagonals of \pentagon_{n+2} that do not involve n+1, and

\bullet a small convex vector \mathsf{B}=(B(P)\,:\,P\text{ is a subpolygon of }\pentagon_{n+2}).

Note that the vector \mathsf{b} of the previous section was concave, and B is convex.

Arkani-Hamed, Bai, He, and Yan [[AHBHY18b](https://arxiv.org/html/2603.03425v3#bib.bibx5)] defined the _associahedron_ 3 3 3 In [[AHBHY18b](https://arxiv.org/html/2603.03425v3#bib.bibx5)] this is denoted \mathcal{A}_{n+2}. We use the subscript n-1 to match the root system A_{n-1} and the dimension of the polytope.\mathrm{ABHY}_{n-1}(\mathsf{A}) to be the (n-1)-dimensional polytope

\displaystyle\{X\in\mathbb{R}^{\operatorname{diag}(\pentagon_{n+2})}\,:\,\displaystyle\,X_{r,s}+X_{r+1,s+1}-X_{r,s+1}-X_{r+1,s}=A_{rs}\text{ for }0\leq r\ll s\leq n
\displaystyle X_{r,s}\geq 0\,\text{ for }\,0\leq r\ll s\leq n+1,\}(9)

where X_{0,n+1}=X_{r,r+1}:=0 for 0\leq r\leq n, and again, r\ll s means that s-r\geq 2.

Arkani-Hamed, Figueiredo, and Vazão [[AHFV25](https://arxiv.org/html/2603.03425v3#bib.bibx8)] defined the following polytope, and conjectured that its face structure matches the combinatorics of Matryoshkas.

###### Definition 4.15.

The _cosmohedron_\mathrm{AFV}_{n-1}(\mathsf{A},\mathsf{B}) is the (n-1)-dimensional polytope

\displaystyle\{X\in\mathbb{R}^{\operatorname{diag}(\pentagon_{n+2})}\,:\,\displaystyle\,X_{r,s}\geq 0\text{ for }0\leq r\ll s\leq n+1,
\displaystyle X_{r,s}+X_{r+1,s+1}-X_{r,s+1}-X_{r+1,s}=A_{rs}\text{ for }0\leq r\ll s\leq n,
\displaystyle\sum_{C\in\mathcal{C}}X_{C}\geq\sum_{P\in\mathcal{S}(\mathcal{C})}B(P)\,:\,\text{ for non-crossing }\mathcal{C}\subset\operatorname{diag}(\pentagon_{n+2})\}(10)

where \mathcal{S}(\mathcal{C}) is the set of polygons in the subdivision of \pentagon_{n+2} induced by \mathcal{C}.

###### Theorem 4.16.

The map \displaystyle X_{r,s}=\sum_{r<i<s}x_{i}-\sum_{r\leq i<<j\leq s}a_{ij} gives linear isomorphisms

1.   1.
between the associahedra \operatorname{Loday}_{n-1}(\mathsf{a}) and \mathrm{ABHY}_{n-1}(\mathsf{A}), and

2.   2.
between the cosmohedra \operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}) and \mathrm{AFV}_{n-1}(\mathsf{A},\mathsf{B}),

where A_{ij}=a_{i,j+1} for each diagonal ij\in\operatorname{diag}(\pentagon_{n+2}) and B(P)=-\varepsilon\,b(P) for sufficiently small \varepsilon>0.

Before proving this, we note that \mathsf{a} has size {\binom{n+1}{2}} whereas \mathsf{A} has size {\binom{n}{2}}; this is because the n entries a_{i-1,i+1} do not enter the definition of \mathrm{AFV}_{n-1}(\mathsf{A},\mathsf{B}). These entries are the coefficients of the points \Delta_{\{i\}}=e_{i} in the Minkowski sum of \operatorname{Loday}_{n-1}(\mathsf{a}), so they only cause a translation of the polytope.

###### Proof of Theorem [4.16](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem16 "Theorem 4.16. ‣ 4.2.4 Isomorphism with the Arkani–Hamed-Figueiredo-Vazão cosmohedron. ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron").

1. A vector x\in\mathbb{R}^{n} determines a vector X in kinematic space that satisfies the linearly independent equations

X_{r,s}+X_{r+1,s+1}-X_{r,s+1}-X_{r+1,s}=a_{r,s+1}=A_{r,s}\qquad\text{ for }0\leq r<<s\leq n.

as one readily verifies. These \binom{n}{2} equations restrict X to a vector space of dimension (n+2)(n-1)/2-\binom{n}{2}=n-1, which is precisely the dimension of the associahedron. The inequalities X_{r-1,s+1}\geq 0 are restatements of the defining inequalities of \operatorname{Loday}_{n-1}(\mathsf{a}) as written in ([4](https://arxiv.org/html/2603.03425v3#S4.E4 "In 4.1 Preliminaries: Associahedra and bracket associahedra ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron")), proving the desired isomorphism. See also [[Ear24](https://arxiv.org/html/2603.03425v3#bib.bibx21)].

2. It remains to translate the inequalities of the cosmohedron \operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}):

\sum_{1\leq i\leq n}d_{\mathcal{C}}(i)\ x_{i}\geq\sum_{ij\in\operatorname{diag}(\pentagon_{n+2})}d_{\mathcal{C}}(ij)\ a_{ij}-\varepsilon\,\sum_{P\in\mathcal{S}}b(P)

for any subdivision \mathcal{S} with chords \mathcal{C}. We rewrite this as

\sum_{C\in\mathcal{C}}\left(\sum_{1\leq i\leq n}d_{\{C\}}(i)\ x_{i}\right)\geq\sum_{C\in\mathcal{C}}\left(\sum_{ij\in\operatorname{diag}(\pentagon_{n+2})}d_{\{C\}}(ij)\ a_{ij}\right)-\varepsilon\,\sum_{P\in\mathcal{S}}b(P).

Now notice that d_{\{C\}}(i) is 1 if C\succ i and 0 otherwise, and d_{\{C\}}(ij) is 1 if C\succeq ij and 0 otherwise, so

X_{C}=\sum_{1\leq i\leq n}d_{\{C\}}(i)\ x_{i}-\sum_{ij\in\operatorname{diag}(\pentagon_{n+2})}d_{\{C\}}(ij)\ a_{ij}.

Thus the inequalities of \operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}) in x-space are equivalent to the inequalities of \mathrm{AFV}_{n-1}(\mathsf{A},\mathsf{B}) in X-space. We note that \mathsf{B} is concave if and only if \mathsf{b} is concave, and \mathsf{B} is made sufficiently small – as [[AHFV25](https://arxiv.org/html/2603.03425v3#bib.bibx8)] requires – by making \varepsilon sufficiently small. ∎

In particular, this proves one of our main results: the correctness of Arkani-Hamed, Figueiredo, and Vazão’s construction.

###### Corollary 4.17.

The face poset of Arkani-Hamed, Figueiredo, and Vazão’s cosmohedron \mathrm{AFV}_{n-1} is anti-isomorphic to the poset of Matryoshkas of an (n+2)-gon.

It also lets us determine precisely what their requirement of a “small enough” B means. Theorem [4.11](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem11 "Theorem 4.11. ‣ 4.2.2 The cosmohedron: proof of correctness ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron") proves that \operatorname{Cosmo}_{n-1}(\mathsf{a},\mathsf{b}) is a cosmohedron for 0<\varepsilon<(\min a_{ij})/24(\max b(P)); it would be interesting to trace the proof carefully and find the optimal constant.

## 5 Face structure and enumeration

The dimension of a face. Recall that M_{\text{not min}} is the number of polygons in a Matryoshka M that are not inclusion-minimal.

###### Lemma 5.1.

The cone of a Matryoshka M in the cosmohedral fan has dimension

\dim\operatorname{cone}(M)=|M_{\text{not min}}|+1.

###### Proof.

If (T,B) is the bracketed tree corresponding to M from Lemma [2.10](https://arxiv.org/html/2603.03425v3#S2.Thmtheorem10 "Lemma 2.10. ‣ 2.3.3 The cosmohedral fan in terms of bracketed trees. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron"), we have \operatorname{cone}(M)=\operatorname{cone}(T,B)=f_{T}(\operatorname{cone}_{T}(B)_{\geq 0}) and its dimension equals the dimension of \operatorname{cone}_{T}(B), which is the number of brackets in B. This equals the number of non-minimal polygons of the Matryoshka M. ∎

Faces of the cosmohedron factor combinatorially, but not geometrically. The faces of the cosmohedron have an interesting factorization property that is subtler and weaker than the analogous property for bracket associahedra.

Bracket associahedra 4 4 4 and, more broadly, generalized permutahedra [[AA23](https://arxiv.org/html/2603.03425v3#bib.bibx1), [Edm70](https://arxiv.org/html/2603.03425v3#bib.bibx22)] have the property that any facet factors canonically as a Cartesian product of two bracket associahedra, as we now describe. Let H be a graph, \beta be a bracket of H, and H/\beta the contraction of \beta in H, whose edges are in bijection with E(H)-E(\beta). Also let b|_{\beta}(\beta^{\prime})=b(\beta^{\prime}) for each bracket \beta^{\prime} of \beta and b/_{\beta}(\beta^{\prime\prime})=b(\beta\cup\beta^{\prime\prime})-b(\beta) for each bracket \beta^{\prime\prime} of H/\beta. Note that if b is concave on \operatorname{Brackets}(H), then b|_{\beta} and b/_{\beta} are concave on \operatorname{Brackets}(\beta) and \operatorname{Brackets}(H/\beta), respectively.

###### Lemma 5.2.

[[AA23](https://arxiv.org/html/2603.03425v3#bib.bibx1)] For any graph H and any bracket \beta of H, the e_{\beta}-maximal facet of \operatorname{BrAssoc}_{H}(\mathsf{b}) is linearly isomorphic to the product of bracket associahedra

(\operatorname{BrAssoc}_{H}(\mathsf{b}))_{\beta}\cong\operatorname{BrAssoc}_{\beta}(\mathsf{b}|_{\beta})\times\operatorname{BrAssoc}_{H/\beta}(\mathsf{b}/_{\beta})\subset\mathbb{R}^{E(\beta)}\times\mathbb{R}^{E(H/\beta)}

Facets, and hence faces, of cosmohedra have an analogous factorization property – but it is significantly weaker. The first statement below was observed in [[AHFV25](https://arxiv.org/html/2603.03425v3#bib.bibx8)], but it holds an unexpected surprise.

###### Proposition 5.3.

Let \mathcal{S} be a subdivision of the (n+2)-gon into polygons P_{1},\ldots,P_{k} where P_{i} is an (n_{i}+2)-gon. Let T(\mathcal{S}) be its dual graph, which is a tree. Then the facet of \operatorname{Cosmo}_{n} corresponding to \mathcal{S} is combinatorially isomorphic to the product

(\operatorname{Cosmo}_{n})_{\mathcal{S}}\cong^{comb}\operatorname{BrAssoc}_{T(\mathcal{S})}\times(\operatorname{Cosmo}_{n_{1}})\times\cdots\times(\operatorname{Cosmo}_{n_{k}}).

However, this factorization does not hold geometrically.

###### Proof.

The face poset of L({\operatorname{Cosmo}_{n}}) is the opposite of the poset M_{n} of Matryoshkas of the (n+2)-gon. The interval [\widehat{0},(\operatorname{Cosmo}_{n})_{\mathcal{S}}] corresponds to the subposet (M_{n})_{\geq{\mathcal{S}}} of Matryoshkas refining \mathcal{S}. Such a Matryoshka is obtained by making k+1 independent choices:

(i) refining each of the k polygons P_{i} internally with its own Matryoshka – the ways of doing this are in order-reversing bijection with the faces of \operatorname{Cosmo}_{n_{i}} – and

(ii) refining \mathcal{S} externally with compatible polygons that are unions of the P_{i}s – the ways of doing this are in order-reversing bijection with the bracketings of T(\mathcal{S}).

These independent choices are compatible with the order relations of the posets, giving the desired result. ∎

The above isomorphism is not, in general, a linear isomorphism of polyhedra. We can already see this in the cosmohedron \operatorname{Cosmo}_{3}. Consider the subdivision of a hexagon by a long and a short diagonal, as illustrated in Figure [18](https://arxiv.org/html/2603.03425v3#S5.F18 "Figure 18 ‣ 5 Face structure and enumeration ‣ Combinatorics of the cosmohedron"). Proposition [5.3](https://arxiv.org/html/2603.03425v3#S5.Thmtheorem3 "Proposition 5.3. ‣ 5 Face structure and enumeration ‣ Combinatorics of the cosmohedron") predicts that this facet is combinatorially isomorphic to the product T(P_{2})\times\operatorname{Cosmo}_{2}\times\operatorname{Cosmo}_{1}\times\operatorname{Cosmo}_{1}=|\times|\times\cdot\times\cdot which is a quadrilateral. We claim that this face is a trapezoid; we can see this in Figure [2](https://arxiv.org/html/2603.03425v3#S1.F2 "Figure 2 ‣ 1.2 Polyhedral properties of the cosmohedron ‣ 1 Introduction ‣ Combinatorics of the cosmohedron"), but let us verify it algebraically. We use Proposition [4.8](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem8 "Proposition 4.8. ‣ 4.2.1 Vertex description ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron") to compute:

\displaystyle c_{M_{1}}=(1,2,6,1)+\varepsilon(0,-3,3,0)\displaystyle c_{N_{1}}=(2,1,6,1)+\varepsilon(-3,0,3,0)
\displaystyle c_{M_{2}}=(2,1,6,1)+\varepsilon(-1,0,3,-2)\displaystyle c_{N_{2}}=(1,2,6,1)+\varepsilon(0,-1,3,-2)

and notice that the opposite sides

c_{M_{1}}-c_{N_{1}}=(1-3\varepsilon)(e_{2}-e_{1}),\qquad c_{M_{2}}-c_{N_{2}}=(1-\varepsilon)(e_{1}-e_{2}),\qquad

are parallel but not equal. Since it is not a parallelogram, this quadrilateral face cannot be expressed as a Cartesian product of two edges.

Figure 18: A trapezoidal face of the cosmohedron.

### 5.1 Enumeration

Having understood the combinatorial structure of the cosmohedron, we are now ready to enumerate its faces. We also compute the face numbers of the _correlatron_ of [[FV26](https://arxiv.org/html/2603.03425v3#bib.bibx23)], which are closely related in an elegant and surprising way.

### 5.2 Counting the faces of the cosmohedron

Our next goal is to compute the number of faces of the cosmohedron in each dimension. We list these numbers in Table [1](https://arxiv.org/html/2603.03425v3#S5.T1 "Table 1 ‣ 5.2 Counting the faces of the cosmohedron ‣ 5 Face structure and enumeration ‣ Combinatorics of the cosmohedron") for n\geq 6.

Table 1: The number of d-codimensional faces of the n th cosmohedron on n+2 points.

It is useful to organize them in the f-polynomial f_{n}(t)=f_{\operatorname{AssocFan}_{n-1}}(t) and the modified f-polynomial g_{n}(t)=[(1+t)f_{n}(t)-1]/t, whose coefficients are the consecutive pair sums of the coefficients of f_{n}(t):

\displaystyle f_{1}(t)=1,\,\,f_{2}(t)=1+2t,\,\,f_{3}(t)=1+10t+10t^{2},\,\,f_{4}(t)=1+44t+114t^{2}+72t^{3},\ldots
\displaystyle g_{1}(t)=1,\,\,g_{2}(t)=3+2t,\,\,g_{3}(t)=11+20t+10t^{2},\,\,g_{4}(t)=45+148t+186t^{2}+72t^{3},\ldots.

###### Theorem 5.4.

The f-polynomials of the cosmohedral fans form the unique sequence \{f_{n}(t)\}_{n\geq 1} of polynomials such that the power series in (\mathbb{R}[t])[x]

\displaystyle x-f_{1}(t)x^{2}-f_{2}(t)x^{3}-f_{3}(t)x^{4}-\ldots
\displaystyle x+g_{1}(t)x^{2}+g_{2}(t)x^{3}+g_{3}(t)x^{3}+\ldots

are compositional inverses, where g_{n}(t)=[(1+t)f_{n}(t)-1]/t. Equivalently, the power series A_{t}(x)=x-f_{1}(t)x^{2}-f_{2}(t)x^{3}+\ldots satisfies

A_{t}\left(\frac{(2t+1)x-(2t+2)x^{2}}{t(1-x)}-\left(1+\frac{1}{t}\right)A_{t}(x)\right)=x.

###### Proof.

Every face F of the cosmohedral fan is contained in a unique minimal face of the associahedral fan A(F). If F corresponds to a Matryoshka M on an (n+2)-gon then A(F) corresponds to the subdivision M_{max} given by the maximal polygons P_{1},\ldots,P_{k} of M. If P_{i} is an (n_{i}+2)-gon then n=n_{1}+\cdots+n_{k}.

Conversely, consider a face A of the associahedral fan, corresponding to a non-trivial subdivision \mathcal{S} of the (n+2)-gon into an (n_{i}+2)-gon P_{i} for 1\leq i\leq k. The faces F of the cosmohedron such that A(F)=A correspond to the Matryoshkas with M_{max}=\mathcal{S}; these are obtained by choosing arbitrary Matryoshkas M_{1},\ldots,M_{k} on P_{1},\ldots,P_{k}, that is, faces F_{1},\ldots,F_{k} on the corresponding cosmohedral fans. Furthermore, we have

\dim F=1+|M_{\text{not min}}|=1+\sum_{i=1}^{k}(|(M_{i})_{\text{not min}}|+1)=1+\sum_{i=1}^{k}\dim F_{i}.

so the f-numbers of the cosmohedron \operatorname{Cosmo}_{n-1} are given by the recurrence

f_{n}(t)=\sum_{F\leq\operatorname{Cosmo}_{n-1}}t^{\operatorname{codim}F}=1+t\sum_{A\lneq\operatorname{Assoc}_{n-1}}f_{n_{1}}(t)\cdots f_{n_{k}}(t).(11)

Adding the full face \operatorname{Assoc}_{n-1} to both sides, we obtain

g_{n}(t)=\sum_{A\leq\operatorname{Assoc}_{n-1}}f_{A}(t)(12)

using the notation of [[AA23](https://arxiv.org/html/2603.03425v3#bib.bibx1)]: if x_{1},x_{2},\ldots, is a sequence and a face A of the associahedron \operatorname{Assoc}_{n-1} is combinatorially a product of associahedra \operatorname{Assoc}_{n_{1}-1}\times\cdots\operatorname{Assoc}_{n_{k}-1}, then we write x_{A}=x_{n_{1}}\cdots x_{n_{k}}.

Arriving at ([12](https://arxiv.org/html/2603.03425v3#S5.E12 "In Proof. ‣ 5.2 Counting the faces of the cosmohedron ‣ 5 Face structure and enumeration ‣ Combinatorics of the cosmohedron")) is a pleasant déjà vu: the associahedral version of Lagrange inversion obtained in [[AA23](https://arxiv.org/html/2603.03425v3#bib.bibx1), Theorem 2.4.3] tells us that this recurrence is precisely equivalent to the assertion that x-f_{1}(t)x^{2}-f_{2}(t)x^{3}-\ldots and x+g_{1}(t)x^{2}+g_{2}(t)x^{3}+\ldots are compositional inverses.

The second formulation follows from the first by a straightforward computation. ∎

A pleasant feature of this proof is that Aguiar and Ardila discovered the associahedral version of Lagrange inversion in [[AA23](https://arxiv.org/html/2603.03425v3#bib.bibx1), Theorem 2.4.3] within the Hopf-algebraic structure of Loday’s associahedra, specifically. Is it coincidental that this is the same associahedron we constructed? Might there be a deeper underlying algebraic explanation? Cosmohedra feature an operadic structure that is the subject of an upcoming paper.

###### Corollary 5.5.

The f-polynomials of the cosmohedral fans are given by the recurrence f_{1}(t)=1 and

f_{n}(t)=1+t\sum_{1^{a_{1}}2^{a_{2}}\cdots\vdash n}\frac{(n+a_{1}+a_{2}+\cdots)!}{(n+1)!}\cdot\frac{f_{1}(t)^{a_{1}}}{a_{1}!}\frac{f_{2}(t)^{a_{2}}}{a_{2}!}\cdots\qquad\text{for }n\geq 2.

summing over the non-trivial partitions of n as a sum n=a_{1}\cdot 1+a_{2}\cdot 2+\cdots+a_{n-1}\cdot(n-1).

###### Proof.

This is a reformulation of the recurrence formula ([11](https://arxiv.org/html/2603.03425v3#S5.E11 "In Proof. ‣ 5.2 Counting the faces of the cosmohedron ‣ 5 Face structure and enumeration ‣ Combinatorics of the cosmohedron")), taking into account that the number of subdivisions of an (n+2)-gon into a_{2} triangles, a_{3} quadrilaterals, etc. is (n+a_{1}+a_{2}+\cdots)!/(n+1)!a_{1}!a_{2}!\cdots[[Sta99](https://arxiv.org/html/2603.03425v3#bib.bibx37), Theorem 5.3.10]. ∎

### 5.3 Asymptotic growth

Facets. The facets of the cosmohedron \operatorname{Cosmo}_{n-1} are in bijection with the non-trivial subdivisions of an (n+2)-gon, or equivalently the non-trivial parenthesizations of a string of n+1 letters. Their number, listed in column 1 of Table [1](https://arxiv.org/html/2603.03425v3#S5.T1 "Table 1 ‣ 5.2 Counting the faces of the cosmohedron ‣ 5 Face structure and enumeration ‣ Combinatorics of the cosmohedron"), is one less than the n-th _Hipparchus–Schröder number_. This is one of the oldest known combinatorial quantities, studied by astronomer and mathematician Hipparchus in the second century B.C.; for a fascinating historical account, see [[Sta99](https://arxiv.org/html/2603.03425v3#bib.bibx37), [Sta97](https://arxiv.org/html/2603.03425v3#bib.bibx36)]. These numbers are sequence A001003 in the Online Encyclopedia of Integer Sequences; their generating function is

x+x^{2}+3x^{3}+11x^{4}+45x^{5}+197x^{6}+\cdots=\frac{1}{4}\left(1+x-\sqrt{1-6x+x^{2}}\right)

The asymptotic growth of the sequence is controlled by 3+\sqrt{8}, the inverse of the smallest root of 1-6x+x^{2}. More precisely [[Knu97](https://arxiv.org/html/2603.03425v3#bib.bibx27)], the number r_{n} of rays grows like

r_{n}\sim\frac{1}{4}\sqrt{(\sqrt{18}-4)/\pi}\cdot(3+\sqrt{8})^{n}\cdot n^{-3/2}.

Vertices. The vertices of the cosmohedron \operatorname{Cosmo}_{n-1} are in bijection with the maximal Matryoshkas on an (n+2)-gon. The following lemma is stated in [[AHFV25](https://arxiv.org/html/2603.03425v3#bib.bibx8)] with a minor typo; we include a proof.

###### Lemma 5.6.

The number m_{n} of maximal Matryoshkas on an (n+2)-gon is given by the recurrence relation

m_{1}=1,\qquad m_{n}=\sum_{k=1}^{n-1}(k+1)m_{k}m_{n-k}\text{ for }n\geq 2.

The formal power series M(x)=\sum_{n\geq 2}m_{n-1}x^{n}=x^{2}+2x^{3}+10x^{4}+72x^{5}+\cdots is given by the differential equation

M^{\prime}(x)=1-\frac{x^{2}}{M(x)}.

###### Proof.

Fix an edge e of the (n+2)-gon. A maximal Matryoshka is obtained by splitting an (n+2)-gon into a (k+2)-gon that contains e and an (n-k+2)-gon, and then putting a maximal Matryoshka in each one of them. There are k+1 possible positions for the outer (k+2)-gon that contains e, and m_{k} and m_{n-k} choices for the two inner Matryoshkas, respectively. The differential equation follows by a straightforward computation. ∎

In particular, this means that M(x) is a D-algebraic series; that is, its derivatives satisfy a polynomial equation with coefficients in \mathbb{C}[x]; in fact, one of the simplest such equations: M(x)=x^{2}+M(x)M^{\prime}(x). Compare this with the algebraic power series for Catalan numbers, which satisfies one of the simplest algebraic equations C(x)=1+xC(x)^{2}.

Melczer [[Mel21](https://arxiv.org/html/2603.03425v3#bib.bibx30)] writes that “_D-algebraic series seem to be on the border between decidability and undecidability, perhaps tipped to the side of undecidability. […] Compared to D-finite functions, not much is known about D-algebraic asymptotics, and in general not much can be said._”. Since M(x) is one of the simplest D-algebraic series, the following conjecture of Kotesovic, based on numerical evidence, is an interesting challenge.

###### Conjecture 5.7.

[[KI14](https://arxiv.org/html/2603.03425v3#bib.bibx26)] The number of Matryoshkas grows asymptotically like

m_{n}\sim c\cdot n!\cdot n^{4}

for c=0.005428\ldots.

All faces. A more ambitious goal would be to find the asymptotics for the total number of faces, or even the limit distribution of the normalized f-vector. We do not know how to do this.

### 5.4 Counting the faces of the correlatron

Another object that naturally arises in physics, and that is intimately related to the cosmohedron, is the _correlatron_ introduced by Figueiredo and Vazão in [[FV26](https://arxiv.org/html/2603.03425v3#bib.bibx23)]. This polytope captures the combinatorics of the correlator, a physical observable that can be computed from the wavefunction. The combinatorial structure of this object is described by the interacting combinatorics of chords (the combinatorics of associahedra) and subpolygons (the combinatorics of cosmohedra).

Figure 19: The three-dimensional correlatron.

The correlatron is an unbounded n-dimensional polyhedron, with faces of two kinds:

1. The bounded faces of the correlatron are labelled by the compatible pairs (\mathcal{C},M), consisting of a set \mathcal{C} of non-crossing chords and a Matryoshka M on an (n+2)-gon. Here \mathcal{C} and M are _compatible_ if no chord in \mathcal{C} intersects a polygon in M; or equivalently, if \mathcal{C} is a subset of the diagonals in the underlying subdivision M_{\text{max}} of M. The bounded face inclusions of the correlatron are given by the reverse inclusion of chords and Matryoshkas:

\{\mathcal{C}^{\prime},M^{\prime}\}\text{ is a face of }\{\mathcal{C},M\}\text{ in the correlatron if and only if }\mathcal{C}\subseteq\mathcal{C}^{\prime}\text{ and }M\leq M^{\prime}.

2. The unbounded faces are labeled by chords \underline{\mathcal{C}}, which we underline to distinguish them easily from the earlier faces. An unbounded face \underline{\mathcal{C}} is contained in \underline{\mathcal{C}^{\prime}} if and only if \mathcal{C}\subseteq\mathcal{C}^{\prime}. An unbounded face \underline{\mathcal{C}} also contains the bounded faces (\mathcal{C}^{\prime},M^{\prime}) such that \mathcal{C}\subseteq\mathcal{C}^{\prime}.

Note in particular that the correlatron has a bottom cosmohedral facet \operatorname{Cosmo}_{n-1} corresponding to the pairs of the form (\emptyset,M). Its unbounded directions look like a cone over \operatorname{Assoc}_{n-1}. Since the unbounded faces of the correlatron are enumerated by the well-known f-vector of the associahedron, we focus on enumerating the bounded part.

Table 2: The number of d-codimensional faces of the n th bounded correlatron.

###### Theorem 5.8.

The reverse f-polynomials of the bounded correlatrons h_{1}(t)=t,h_{2}(t)=3t+4t^{2},h_{3}(t)=11t+35t^{2}+25t^{3},\ldots form the unique sequence \{h_{n}(t)\}_{n\geq 1} of polynonmials such that the power series in (\mathbb{R}[t])[x]

\displaystyle x+h_{1}(t)x^{2}+h_{2}(t)x^{3}+h_{3}(t)x^{4}+\ldots
\displaystyle x-tg_{1}(t)x^{2}-tg_{2}(t)x^{3}-tg_{3}(t)x^{3}+\ldots

are compositional inverses, where g_{n}(t)=[(1+t)f_{n}(t)-1]/t and f_{n} is the f-polynomial of the n th cosmohedral fan. Equivalently, the power series H_{t}(x)=x+h_{1}(t)x^{2}+h_{2}(t)x^{3}+\ldots satisfies

H_{t}\left(\frac{(t+1)x^{2}-tx}{1-x}+(t+1)A_{t}(x)\right)=x

where A_{t}(x) is the generating function for the f-vectors of cosmohedra in Theorem [5.4](https://arxiv.org/html/2603.03425v3#S5.Thmtheorem4 "Theorem 5.4. ‣ 5.2 Counting the faces of the cosmohedron ‣ 5 Face structure and enumeration ‣ Combinatorics of the cosmohedron").

###### Proof.

Every face F of the truncated correlatron corresponds to a pair (\mathcal{C},M) where M is non-trivial. Let \mathcal{C} create a subdivision into polygons P_{1},\ldots,P_{k}, and let A(F) be the corresponding face of the associahedal fan corresponding to the set of chords \mathcal{C}. A Matryoshka M compatible with \mathcal{C} restricts to Matryoshkas M_{1},\ldots,M_{k} on P_{1},\ldots,P_{k}. One can almost recover M from the M_{i}s: the only additional information is whether M contains the polygon P_{i} (in which case we let \varepsilon_{i}=1) or it does not (in which case \varepsilon_{i}=0). We note that if M_{i} is the trivial Matryoshka on P_{i}, then \varepsilon_{i}=1 necessarily; if M_{i} is non-trivial, then \varepsilon_{i} can be chosen freely. Finally notice that

\operatorname{codim}F=|\mathcal{C}|+|M_{\text{non-min}}|=k+\sum_{i=1}^{k}|(P_{i})_{\text{non-min}}|+\sum_{i=1}^{k}(\varepsilon_{i}-1)

Therefore

\displaystyle h_{n}(t)\displaystyle=\displaystyle\sum_{F\leq\operatorname{BdCorrela}_{n}}t^{\operatorname{codim}F}
\displaystyle=\displaystyle\sum_{\begin{subarray}{c}A\leq\operatorname{Assoc}_{n}\\
\mathcal{S}(A)=\{P_{1},\ldots,P_{k}\}\end{subarray}}t^{k}\sum_{\begin{subarray}{c}M_{1},\ldots,M_{k}\\
\varepsilon_{1},\ldots,\varepsilon_{k}\end{subarray}}\left(\prod_{i=1}^{k}t^{|(P_{i})_{\text{non-min}}|+\varepsilon_{i}-1}\right)
\displaystyle=\displaystyle\sum_{A\leq\operatorname{Assoc}_{n}}t^{|\mathcal{C}|}\prod_{i=1}^{k}\left(1+\sum_{M_{i}\neq\{P_{i}\}}t^{|(P_{i})_{\text{non-min}}|}\left(1+\frac{1}{t}\right)\right)
\displaystyle=\displaystyle\sum_{A\leq\operatorname{Assoc}_{n}}t^{|\mathcal{C}|}\prod_{i=1}^{k}\left(1+(f_{P_{i}}(t)-1)\left(1+\frac{1}{t}\right)\right)
\displaystyle=\displaystyle\sum_{A\leq\operatorname{Assoc}_{n}}t\,\prod_{i=1}^{k}\left((1+t)f_{P_{i}}(t)-1\right)
\displaystyle=\displaystyle\sum_{A\leq\operatorname{Assoc}_{n}}tg_{A}(t),

again using the notation of [[AA23](https://arxiv.org/html/2603.03425v3#bib.bibx1)]. The associahedral version of Lagrange inversion [[AA23](https://arxiv.org/html/2603.03425v3#bib.bibx1), Theorem 2.4.3] then tells us that x-tg_{1}(t)x^{2}-tg_{2}(t)x^{3}-tg_{3}(t)x^{4}+\ldots and x+h_{1}(t)x^{2}+h_{2}(t)x^{3}+h_{3}(t)x^{3}+\ldots are compositional inverses.

The second formulation follows from a straightforward computation. ∎

Note the remarkable similarity between the formulas in Theorems [5.4](https://arxiv.org/html/2603.03425v3#S5.Thmtheorem4 "Theorem 5.4. ‣ 5.2 Counting the faces of the cosmohedron ‣ 5 Face structure and enumeration ‣ Combinatorics of the cosmohedron") and [5.8](https://arxiv.org/html/2603.03425v3#S5.Thmtheorem8 "Theorem 5.8. ‣ 5.4 Counting the faces of the correlatron ‣ 5 Face structure and enumeration ‣ Combinatorics of the cosmohedron") for the f-vectors of the cosmohedron and bounded correlatron. Is this an indication of something deeper?

## 6 Elusive qualities of the cosmohedron

The cosmohedron is related to several beautiful polytopes in the literature, either through direct connections or indirect analogies. These include the permutahedron [[BS13](https://arxiv.org/html/2603.03425v3#bib.bibx14), [Sch13](https://arxiv.org/html/2603.03425v3#bib.bibx34)], the associahedron [[Lod04](https://arxiv.org/html/2603.03425v3#bib.bibx29), [Sta63](https://arxiv.org/html/2603.03425v3#bib.bibx35)], graph associahedra [[CD06](https://arxiv.org/html/2603.03425v3#bib.bibx15)], the nested permutahedron [[CL22](https://arxiv.org/html/2603.03425v3#bib.bibx16)], the permutoassociahedron [[CL23](https://arxiv.org/html/2603.03425v3#bib.bibx17)], and the permutonestohedron [[Gai15](https://arxiv.org/html/2603.03425v3#bib.bibx24)], among others.

However, there are numerous ways in which the combinatorial structure of the cosmohedron is more subtle than its predecessors. In turn, this makes it significantly harder to prove the correctness of this construction. In this section, we discuss some of the intricacies that arise – and that remain in future polytopes of interest.

The cosmohedron is not very symmetric. The permutahedron, the nested permutahedron, the permuto-associahedron, and the permutonestohedron are symmetric under the natural action of S_{n} onto \mathbb{R}^{n}. This allows us to express them as the convex hull of the S_{n}-orbit of a point, or of a carefully placed permutahedron or associahedron, respectively. The cosmohedron only has a natural action of the dihedral group.

The cosmohedron is made of an associahedron’s worth of different bracket associahedra. It requires much greater care to choose a suitable realization of each bracket associahedron that are compatible with each other, and then place each one of them at the correct position, in order to make sure that they glue correctly into a well-behaved global polytope. This is a subtle matter.

The cosmohedron is not simplicial. This is not a surprise; none of the polytopes on the above list are simplicial. Proposition [5.3](https://arxiv.org/html/2603.03425v3#S5.Thmtheorem3 "Proposition 5.3. ‣ 5 Face structure and enumeration ‣ Combinatorics of the cosmohedron") shows that the only simplicial faces of the cosmohedron are the vertices and edges.

The cosmohedron is not simple. This is more surprising, because all the other polytopes in the list above are simple. The only simple vertices of the cosmohedron are those whose Matryoshka is nested linearly. To see this, let M=(T,B) be a maximal Matryoshka. The facets of \operatorname{cone}(M) correspond to the facets of \operatorname{cone}_{T}(B)_{\geq 0}. These in turn correspond to the n-2 facets of \operatorname{cone}_{T}(B) (which is simplicial), and the facets of the form y_{e}\geq 0 for a minimal bracket in B, which is an edge in B, and corresponds to a quadrilateral in M. Therefore

(\text{number of facets of }\operatorname{cone}(M))=n-2+(\text{number of quadrilaterals in }M).

This can be any number between n-1 and \lfloor 3n/2\rfloor-2. It equals n-1 precisely if M has exactly one quadrilateral, which happens precisely when M is linearly nested.

Why is this a challenge? The fact that polytopes like (bracket) associahedra are simple allows for a flexible construction. One can start with an easy polytope (in these cases a simplex), and then create new facets one at a time by making cuts that are “not too deep”, so that the facets meet transversally. The required care only involves depth inequalities; see [[Lod04](https://arxiv.org/html/2603.03425v3#bib.bibx29), [Dev09](https://arxiv.org/html/2603.03425v3#bib.bibx20), [PPP23](https://arxiv.org/html/2603.03425v3#bib.bibx33)].

We do construct the cosmohedron through a chiseling of faces of the associahedron, but we have to be much more careful now: the chiseling has to line up perfectly to create the required non-transversal crossings. The required care now involves many equalities as well as inequalities. In \mathrm{AFV}_{n-1}, these equalities manifest in the fact that the linear forms \sum_{C\in\mathcal{C}}X_{C} satisfy many linear relations between them. In \operatorname{Cosmo}_{n-1}, they manifest in the many requirements b_{B}=b_{B^{\prime}} in Remark [4.12](https://arxiv.org/html/2603.03425v3#S4.Thmtheorem12 "Remark 4.12. ‣ 4.2.2 The cosmohedron: proof of correctness ‣ 4.2 The cosmohedron ‣ 4 The cosmohedron ‣ Combinatorics of the cosmohedron").

The cosmohedron depends delicately on the chiseling. To further illustrate the point above, recall that in our construction of the cosmohedron, we chiseled a Devadoss bracket associahedron at each vertex of the associahedron. This choice is important. For instance, despite their usefulness in other settings, the Postnikov bracket associahedra do not produce the cosmohedron upon chiseling.

We don’t know the cosmohedron to be a deformation of an easy polytope. The _deformations of the permutahedron_, also known as _generalized permutahedra_ or _polymatroids_, are obtained from the permutahedron by moving facets without pushing them past vertices. They are characterized by having the same edge directions as the permutahedron, namely e_{i}-e_{j} for i\neq j[[ACEP20](https://arxiv.org/html/2603.03425v3#bib.bibx2), Prop. 2.6], [[Pos09](https://arxiv.org/html/2603.03425v3#bib.bibx32)].

Generalizing this setup, Castillo and Liu [[CL22](https://arxiv.org/html/2603.03425v3#bib.bibx16)] defined the _nested permutahedron_\operatorname{Perm}(\alpha,\beta) by placing an (n-1)-permutahedron at each corner of the n-permutahedron. They studied its deformations, and showed that Kapranov’s _permutoassociahedron_[[Kap93](https://arxiv.org/html/2603.03425v3#bib.bibx25)] – which has an (n-1)-associahedron at each vertex of the n-permutahedron, can be obtained as a deformation of \operatorname{Perm}(\alpha,\beta).

The cosmohedron would seem to fit in the same framework; our proof of Theorem [3.6](https://arxiv.org/html/2603.03425v3#S3.Thmtheorem6 "Theorem 3.6. ‣ 3.3 The cosmohedral fan is a complete fan ‣ 3 Proof of correctness of the cosmohedral fan ‣ Combinatorics of the cosmohedron") even shows that the cosmohedron has the same edge directions as \operatorname{Perm}(\alpha,\beta), namely e_{i}-e_{j} and e_{i}-e_{j}+e_{k}-e_{l} for distinct i,j,k,l. This makes it natural to wonder: Is our cosmohedron a deformation of the nested permutahedron? It is not.

Dually, we claim that our cosmohedral fan is not a coarsening of the nested braid fan \operatorname{Br}_{n}^{2} of Castillo and Liu [[CL22](https://arxiv.org/html/2603.03425v3#bib.bibx16)]. The n!\cdot(n-1)! facets of \operatorname{Br}_{n}^{2} correspond to the choices of an order \sigma\in S_{n} of the coordinates of x\in\mathbb{R}^{n} and an order \tau\in S_{n-1} of their consecutive differences. For example, the set of points x\in\mathbb{R}^{4} with

x_{3}>x_{1}>x_{4}>x_{2} and x_{3}-x_{1}>x_{1}-x_{4}>x_{4}-x_{2}>0

is a maximal open cone of \operatorname{Br}_{4}^{2}. Such a point must satisfy two inequalities x_{3}-x_{4}>x_{3}-x_{1} and x_{1}-x_{2}>0 of the Matryoshkal \operatorname{cone}(M_{2}) of Figure [9](https://arxiv.org/html/2603.03425v3#S2.F9 "Figure 9 ‣ 2.3.3 The cosmohedral fan in terms of bracketed trees. ‣ 2.3 Bracketed trees and the cosmohedral fan ‣ 2 Matryoshkas, bracketed trees, and the cosmohedral fan ‣ Combinatorics of the cosmohedron"), but it may or may not satisfy the third inequality x_{3}-x_{1}>x_{1}-x_{2}. Therefore, this maximal cone of \operatorname{Br}_{4}^{2} intersects \operatorname{cone}(M_{2}) full-dimensionally, but is not contained in it. This proves our claim.

Why is this a challenge? Generalized permutahedra are in bijection with submodular functions on the subsets of [n]. Thus to construct them, it is sufficient to choose the correct submodular function, and to use the theory of polymatroids and generalized permutahedra to prove the correct combinatorial structure; see [[AA23](https://arxiv.org/html/2603.03425v3#bib.bibx1), [Pos09](https://arxiv.org/html/2603.03425v3#bib.bibx32)].

We do not know an analog to the permutahedron for this setting. It would be very interesting to find a nice polytope such that (i) it can be deformed to obtain the cosmohedron, (ii) its deformation cone can be nicely described combinatorially, and (iii) the face structure of its deformations is combinatorially well-behaved.

We don’t know the cosmohedron to be a Minkowski sum of easy polytopes. The permutahedron, Loday’s associahedron, and Postnikov’s graph and bracket associahedra can all be expressed as Minkowski sums of faces of the standard simplex. This gives a simple proof that their face structure– the common refinement of the easy face structures of the summands – is correct. For details, see [[Pos09](https://arxiv.org/html/2603.03425v3#bib.bibx32)].

We do not know how to do this in our setting. It would be interesting to find a realization of the cosmohedron that can be decomposed easily and usefully into a Minkowski sum.

The faces of the cosmohedron only factor combinatorially. As discussed in Lemma [5.2](https://arxiv.org/html/2603.03425v3#S5.Thmtheorem2 "Lemma 5.2. ‣ 5 Face structure and enumeration ‣ Combinatorics of the cosmohedron"), faces of a generalized permutahedra factor into smaller generalized permutahedra [[AA23](https://arxiv.org/html/2603.03425v3#bib.bibx1), [Edm70](https://arxiv.org/html/2603.03425v3#bib.bibx22)]. For permutahedra, associahedra, and graph associahedra, any facet of a polytope in these families is a product of smaller polytopes in the same family. This gives an inductive strategy to prove their combinatorial structure, introduced for associahedra in [[MSS02](https://arxiv.org/html/2603.03425v3#bib.bibx31)]. One simply needs to guarantee that the facets are placed correctly; again, because these polytopes are simple, this step only requires verifying the appropriate depth inequalities. Once the facets are laid out correctly, thanks to their factorization, the rest of the face structure can be proved recursively.

In our realization of the cosmohedron, each one of the facets is combinatorially, but not geometrically a product of cosmohedra and a bracket associahedron. Therefore, its combinatorial structure cannot be proved using this inductive strategy. It would be interesting to find a realization of the cosmohedron where the face factorization of Proposition [5.3](https://arxiv.org/html/2603.03425v3#S5.Thmtheorem3 "Proposition 5.3. ‣ 5 Face structure and enumeration ‣ Combinatorics of the cosmohedron") holds geometrically.

Outlook: chiseled polytopes. In an upcoming paper, we present a general family of _\mathcal{X} in Y_ polytopes that contains all of the polytopes above, and numerous other polytopes of interest, old and new. Our proof strategy for the combinatorial structure of the cosmohedron is very well suited to this general setting. In the Appendix, we provide some physical motivation for this general construction, by sketching a new motivating example: we introduce the _one-loop \mathcal{U}-polytope_, which captures the combinatorics of all one-loop cubic (trivalent) graphs, as well as the Feynman polytope for the Symanzik \mathcal{U}-polynomial of each graph, all in one polytope.

## 7 Acknowledgements

The work of Federico Ardila–Mantilla was supported by the National Science Foundation (Grant DMS 2154279) in the United States, a Wolfson Fellowship of The Royal Society in the United Kingdom, and the Friends of the Institute for Advanced Study Fund at the IAS. The work of Nima Arkani-Hamed was supported by the DOE (Grant No. DE-SC0009988), the Simons Collaboration on Celestial Holography, the European Union (ERC, UNIVERSE PLUS, 101118787), and the Carl B. Feinberg cross-disciplinary program in innovation at the IAS. The work of Carolina Figueiredo was supported by FCT/Portugal (Grant No. 2023.01221.BD). The work of Francisco Vazão was supported by the Jonathan M. Nelson Center for Collaborative Research. Views and opinions expressed are those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them.

## Appendix A Appendix: \mathcal{X} in Y polytopes

The main focus of this paper has been to explore and prove the mathematical properties of cosmohedra, which encode the cosmological wavefunction of \operatorname{Tr}(\phi^{3}) theory. We use this appendix to discuss how the principles that guided the construction of cosmohedra are of much broader relevance.

### A.1 The cosmohedron as a guiding example for \mathcal{X} in Y polytopes

As we have learned by now, cosmohedra encode geometrically the combinatorics of maximal Matryoshkas of an n-gon. One way to organize all such Matryoshkas is by dividing them into classes labeled by the underlying collection of triangles – i.e. the underlying triangulation. This is what naturally occurs in physics. To write down the cosmological wavefunction, we start by listing all possible cubic (trivalent) graphs, which correspond to all possible triangulations. Then, for each cubic graph, we consider all possible ways of finding maximal collections of nested subgraphs, which give the possible maximal Matryoshkas compatible with the respective underlying triangulation. The cosmohedron is then the object that puts both combinatorial problems together. In particular, to find the cosmohedron we first needed to find two other polytopes:

1.   Y.
The underlying polytope whose vertices label all possible triangulations/diagrams – which at tree-level is simply the associahedron;

2.   \mathcal{X}.
A class of polytopes, one associated to each triangulation, whose vertices label all possible Matryoshkas compatible with the respective triangulation – which is done by the bracket associahedon X_{T} for each triangulation T.

Once we have these two objects under control, we “blow up” each vertex v_{T} of Y into the respective \mathcal{X} polytope X_{T}, via the chiseling procedure described earlier. At the level of the normal fans, we subdivide each cone of the fan of Y using the corresponding normal fan in the family \mathcal{X}.

Therefore, for any given combinatorial problem, where on top of listing all possible diagrams (Y), we also want to keep track of some further information X_{G} of each diagram G (\mathcal{X}). Once we find the polytopes for \mathcal{X} and Y separately, we can seek to define the polytope that captures the full combinatorics (\mathcal{X} in Y) by the procedure described above.

As it turns out, such situations are ubiquitous in physics, and not just in the mathematical settings discussed in Section [6](https://arxiv.org/html/2603.03425v3#S6 "6 Elusive qualities of the cosmohedron ‣ Combinatorics of the cosmohedron"). In this appendix, we will discuss a particular example, simply to illustrate how this construction generalizes to other contexts.

### A.2 The \mathcal{U}-polytope for one-loop-integrated amplitudes

At loop-level, the amplitudes are expressed as a sum of _Feynman integrals_, which often diverge at the boundaries of the integration domain. Thus, it is important to record the asymptotic behavior of the integrands. The most transparent way to do this is by using the projective representation of Feynman integrals in terms of _Schwinger parameters_\alpha_{e} – one for each edge of the Feynman diagram G – where they take the following form:

\mathcal{I}_{G}=\Gamma(d)\int_{\mathbb{R}\mathbb{P}_{>0}^{E-1}}\frac{\prod_{e\in G}d\alpha_{e}}{\text{GL}(1)}\,\frac{1}{\mathcal{U}_{G}^{D/2-d}\mathcal{F}_{G}^{d}}.

Here d=E-LD/2 is the _superficial degree of divergence_ of a graph G with E edges and genus L (number of independent loops) in spacetime dimension D, and \mathcal{U}_{G} and \mathcal{F}_{G} are the first and second _Symanzik polynomials_ of G:

\mathcal{U}_{G}=\sum_{T\text{ tree of }G}\,\prod_{e\notin T}\alpha_{e},\qquad\mathcal{F}_{G}=\sum_{T\text{ 2-tree of }G}\,\prod_{e\notin T}\alpha_{e}.

The sum in \mathcal{U}_{G} is over all spanning trees of G, and the sum in \mathcal{F}_{G} is over all spanning _2-trees_ T_{1}\sqcup T_{2}, where T_{1} and T_{2} are two vertex-disjoint trees that span all vertices of G. To correctly describe the divergences of \mathcal{I}_{G}, one needs to characterize the asymptotic behavior of \mathcal{U}_{G} and \mathcal{F}_{G}.

Therefore, at one-loop level, we are looking for an object that captures not only all possible cubic (trivalent) one-loop graphs G (collected in a polytope Y), but also the possible leading monomials in \mathcal{U}_{G} and \mathcal{F}_{G} for each graph G (collected in a family of polytopes \mathcal{X}). As it turns out, the polytopes Y and \mathcal{X} polytopes are known in the physics and mathematics literature [[AH](https://arxiv.org/html/2603.03425v3#bib.bibx3), [AHHM22](https://arxiv.org/html/2603.03425v3#bib.bibx9), [BR25](https://arxiv.org/html/2603.03425v3#bib.bibx13), [Edm70](https://arxiv.org/html/2603.03425v3#bib.bibx22)]. We will show how to combine them together into a single chiseled polytope. This answers a question posed in [[AHFS+25](https://arxiv.org/html/2603.03425v3#bib.bibx7)].

In this Appendix, we illustrate the essential ideas behind the construction in a simple but already interesting case: the asymptotic behaviors of the \mathcal{U}-polynomials in the one-loop setting. In an upcoming paper, we will treat the general setting of higher loops, and consider both Symanzik polynomials \mathcal{F} and \mathcal{U}.

#### A.2.1 The Y polytope: the one-loop associahedron

The _one-loop associahedron_ is a polytope designed to parameterize the triangulations of a polygon \hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n} with one puncture. Consider an n-gon with vertices cyclically labeled 1,2,\ldots,n and a puncture labeled O. We consider curves of two kinds:

1.   1.
_chords_ x_{i,j} that go from i to j while keeping the puncture O to the right, for 1\leq i,j\leq n. (Note that we have _diagonals_ x_{i,j}\neq x_{j.i} for each i\neq j, and _tadpoles_ x_{i,i} for each i.)

2.   2.
_spokes_ y_{k} and \tilde{y}_{k} that spiral from k into the puncture clockwise and counterclockwise, respectively. (For simplicity, we draw straight lines from k to O tagged in blue and red, respectively.)

We are interested in the triangulations of \hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n} using these curves. Notice that every triangulation uses exactly n curves, and contains at least one spoke y_{k} or \tilde{y}_{k}. Also, notice that a clockwise and a counterclockwise spiral always intersect, so every triangulation is monochromatic. The left panels of Figures [20](https://arxiv.org/html/2603.03425v3#A1.F20 "Figure 20 ‣ A.2.3 Towards a general construction: 𝒳 in 𝑌 fans in dimensions 2 and 3 ‣ A.2 The 𝒰-polytope for one-loop-integrated amplitudes ‣ Appendix A Appendix: 𝒳 in Y polytopes ‣ Combinatorics of the cosmohedron") and [21](https://arxiv.org/html/2603.03425v3#A1.F21 "Figure 21 ‣ Example A.2. ‣ A.2.3 Towards a general construction: 𝒳 in 𝑌 fans in dimensions 2 and 3 ‣ A.2 The 𝒰-polytope for one-loop-integrated amplitudes ‣ Appendix A Appendix: 𝒳 in Y polytopes ‣ Combinatorics of the cosmohedron") show the 6 triangulations (and six coarser subdivisions) for n=2, and the 10 blue triangulations for n=3, respectively

The _loop associahedron_ was first introduced in [[AH](https://arxiv.org/html/2603.03425v3#bib.bibx3)]; it is a variant of the cluster polytopes in [[AHHST22](https://arxiv.org/html/2603.03425v3#bib.bibx10)] that includes the _tadpole_ curves, x_{i,i}. For the punctured n-gon, we get an n-dimensional simplicial polytope whose vertices correspond to the triangulations of \hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n}. Its faces correspond to the partial triangulations or subdivisions of \hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n}, and the containment relations of the faces match the reverse inclusion relations of the sets of curves. The left panels of Figures [22](https://arxiv.org/html/2603.03425v3#A1.F22 "Figure 22 ‣ A.2.5 The 𝒰-polytope ‣ A.2 The 𝒰-polytope for one-loop-integrated amplitudes ‣ Appendix A Appendix: 𝒳 in Y polytopes ‣ Combinatorics of the cosmohedron") and [23](https://arxiv.org/html/2603.03425v3#A1.F23 "Figure 23 ‣ Example A.3. ‣ A.2.5 The 𝒰-polytope ‣ A.2 The 𝒰-polytope for one-loop-integrated amplitudes ‣ Appendix A Appendix: 𝒳 in Y polytopes ‣ Combinatorics of the cosmohedron") respectively illustrate the one-loop associahedra for \hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{2} and \hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{3}.

#### A.2.2 The \mathcal{X} polytopes: Feynman polytopes

The asymptotic behavior of the Symanzik polynomials \mathcal{U}_{G} and \mathcal{F}_{G} is fully characterized by their Newton polytopes in \mathbb{R}^{E(G)}:

U_{G}=\operatorname{Newt}(\mathcal{U}_{G}),\qquad F_{G}=\operatorname{Newt}(\mathcal{F}_{G})

called _Feynman polytopes_ and studied in [[AHHM22](https://arxiv.org/html/2603.03425v3#bib.bibx9)].

Since \mathcal{U}_{G} and \mathcal{F}_{G} are squarefree, U_{G} and F_{G} are 0-1 polytopes. The exponents of \mathcal{U}_{G} (resp. \mathcal{F}_{G}) are the complements of the spanning trees (resp. 2-trees) of graph G, which are the bases of the matroid M(G) (resp. the truncation matroid \operatorname{Tr}M(G)). It follows that

U_{G}=P_{M(G)^{\perp}},\qquad F_{G}=P_{(\operatorname{Tr}M(G))^{\perp}}

where P_{M} denotes the matroid (basis) polytope of M and ⟂ denotes matroid duality.

In this appendix, for clarity and simplicity, we focus on the first Symanzik polynomials in the one-loop setting. When G only has one cycle, the polynomial \mathcal{U}_{G} is simply the sum of all variables \alpha_{e} associated to the edges e of that cycle, and the Feynman polytope is just the standard simplex on those edges. We summarize this by writing that U_{G}=\Delta_{\operatorname{Cycle}(G)} in the one-loop setting. Note that, unlike in the setting of cosmohedra, this polytope is usually low-dimensional.

#### A.2.3 Towards a general construction: \mathcal{X} in Y fans in dimensions 2 and 3

To motivate the general construction, let us first carry it out in detail in dimensions 2 and 3, beginning at the level of fans. The fan of the one-loop associahedron has a ray for each curve of \hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n}, and a face \operatorname{cone}(\mathcal{C}) for each collection of \mathcal{C}\subset\operatorname{Curves}(\hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n}) of compatible curves, or equivalently, each subdivision of \hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n}. We need to subdivide the cones of this fan to keep track of the dominant terms of the corresponding \mathcal{U}-polynomials.

Figure 20:  (a) Fan of the two-dimensional one-loop associahedron. (b) Refinement of a cone in the fan to encode the dominating terms of the U-polynomial. 

###### Example A.1.

For n=2 there are three blue and three red triangulations of \hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{2}. Their cones naturally fit into a complete fan in the plane shown in the left panel of Figure [20](https://arxiv.org/html/2603.03425v3#A1.F20 "Figure 20 ‣ A.2.3 Towards a general construction: 𝒳 in 𝑌 fans in dimensions 2 and 3 ‣ A.2 The 𝒰-polytope for one-loop-integrated amplitudes ‣ Appendix A Appendix: 𝒳 in Y polytopes ‣ Combinatorics of the cosmohedron"). The only diagrams with a non-trivial \mathcal{U} polynomial are the ones with two spokes, corresponding to the positive and negative quadrants of the fan. The dual diagram of these triangulations is the bubble diagram with two parallel edges \alpha_{1} and \alpha_{2}; its \mathcal{U}-polynomial is simply: \mathcal{U}\left(\begin{gathered}\hbox to46.12pt{\vbox to23.36pt{\pgfpicture\makeatletter\hbox{\hskip 23.06151pt\lower-11.68092pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}\lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.6pt} \lx@inpgf@ignorespaces
{}{{\lx@inpgf@ignorespaces}}\lx@inpgf@ignorespaces\hbox{\hbox{{\lxSVG@begingroup@{_scopebegin=1} {{}{{}}{}{}}
\lxSVG@closescope }}}
{}{{\lx@inpgf@ignorespaces}}\lx@inpgf@ignorespaces\hbox{\hbox{{\lxSVG@begingroup@{_scopebegin=1} {{}{{}}{}{}}
\lxSVG@closescope }}}
{}{{\lx@inpgf@ignorespaces}}\lx@inpgf@ignorespaces\hbox{\hbox{{\lxSVG@begingroup@{_scopebegin=1} {{}{{}}{}{}}
\lxSVG@closescope }}}
{}{{\lx@inpgf@ignorespaces}}\lx@inpgf@ignorespaces\hbox{\hbox{{\lxSVG@begingroup@{_scopebegin=1} {{}{{}}{}{}}
\lxSVG@closescope }}}
{{}}{}{{}}{}{{}}
{}{}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M -31.5 0 L -15.75 0}{fill:none} \lx@inpgf@ignorespaces
{{}}{}{{}}{}{{}}
{}{}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 31.5 0 L 15.75 0}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 0 M 15.75 0 C 15.75 8.7 8.7 15.75 0 15.75 C -8.7 15.75 -15.75 8.7 -15.75 0 C -15.75 -8.7 -8.7 -15.75 0 -15.75 C 8.7 -15.75 15.75 -8.7 15.75 0 Z M 0 0}{fill:none} \lx@inpgf@ignorespaces
{{}}{{\lx@inpgf@ignorespaces}}\lx@inpgf@ignorespaces\hbox{\hbox{{\lxSVG@begingroup@{_scopebegin=1} {{}{}{{
{}{}}}{
{}{}}
{{}{{\lx@inpgf@ignorespaces}}}{{}{\lx@inpgf@ignorespaces}}{}{{}{\lx@inpgf@ignorespaces}}{}{}{}{}{}
{\lx@inpgf@ignorespaces
}{{{{\lx@inpgf@ignorespaces}}\lxSVG@begingroup@{_scopebegin=1} \lxSVG@transformcm{0.7}{0.0}{0.0}{0.7}{-3.80908pt}{5.84668pt}\lxSVG@begingroup@{transform=matrix(0.7 0.0 0.0 0.7 -5.27 8.09)} \pgfsys@hbox{61}\lxSVG@closescope }}}
\lxSVG@closescope }}}
{{}}{{\lx@inpgf@ignorespaces}}\lx@inpgf@ignorespaces\hbox{\hbox{{\lxSVG@begingroup@{_scopebegin=1} {{}{}{{
{}{}}}{
{}{}}
{{}{{\lx@inpgf@ignorespaces}}}{{}{\lx@inpgf@ignorespaces}}{}{{}{\lx@inpgf@ignorespaces}}{}{}{}{}{}
{\lx@inpgf@ignorespaces
}{{{{\lx@inpgf@ignorespaces}}\lxSVG@begingroup@{_scopebegin=1} \lxSVG@transformcm{0.7}{0.0}{0.0}{0.7}{-3.80908pt}{-7.81055pt}\lxSVG@begingroup@{transform=matrix(0.7 0.0 0.0 0.7 -5.27 -10.81)} \pgfsys@hbox{61}\lxSVG@closescope }}}
\lxSVG@closescope }}}
\lxSVG@closescope {{
{}{}{}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}}{{{\lx@inpgf@ignorespaces}}{{}}}}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}\end{gathered}\right)=\alpha_{1}+\alpha_{2},.

We wish to divide each such cone into two subcones, depending on which term of \mathcal{U} dominates: \alpha_{1} or \alpha_{2}. We can easily do this by subdividing them in half, adding the rays g_{1,2}=g_{1}+g_{2} and \tilde{g}_{1,2}=\tilde{g}_{1}+\tilde{g}_{2} respectively. Here g_{i} are the usual _g-vectors_ for the spokes curves, which are the rays of the loop associahedral fan. This subdivision is shown in the right panel of Figure [20](https://arxiv.org/html/2603.03425v3#A1.F20 "Figure 20 ‣ A.2.3 Towards a general construction: 𝒳 in 𝑌 fans in dimensions 2 and 3 ‣ A.2 The 𝒰-polytope for one-loop-integrated amplitudes ‣ Appendix A Appendix: 𝒳 in Y polytopes ‣ Combinatorics of the cosmohedron").

###### Example A.2.

For n=3 there are 10 blue (and 10 red) triangulations of the punctured disk with three marked points on the boundary \hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{3}. They label the maximal cones of the fan shown on the left of Figure [21](https://arxiv.org/html/2603.03425v3#A1.F21 "Figure 21 ‣ Example A.2. ‣ A.2.3 Towards a general construction: 𝒳 in 𝑌 fans in dimensions 2 and 3 ‣ A.2 The 𝒰-polytope for one-loop-integrated amplitudes ‣ Appendix A Appendix: 𝒳 in Y polytopes ‣ Combinatorics of the cosmohedron"); this is the blue half of the one-loop associahedral fan. The red half is isomorphic, and is glued to this one along the boundary in the opposite orientation. Together they form a complete three-dimensional fan.

![Image 4: Refer to caption](https://arxiv.org/html/2603.03425v3/New1Loop3ptFan.png)

Figure 21: (a) Half-fan of the three-dimensional one-loop associahedron. (b) Half-fan of the three-dimensional one-loop \mathcal{U}-polytope. 

Again, for each triangulation \mathcal{C} with two spokes y_{i} and y_{j}, we introduce a new ray g_{ij}=g_{i}+g_{j} that subdivides the corresponding cone in half. But now there is also a triangulation with three spokes y_{1},y_{2},y_{3}, with \mathcal{U} polynomial \alpha_{1}+\alpha_{2}+\alpha_{3}. We subdivide its face \operatorname{cone}(g_{1},g_{2},g_{3}) by introducing four new rays g_{1,2}=g_{1}+g_{2}, g_{2,3}=g_{2}+g_{3},g_{1,3}=g_{1}+g_{3}, and that span three full-dimensional subcones: \operatorname{cone}_{1}=\operatorname{cone}(g_{1},g_{1,2},g_{1,3},g_{1,2,3}), \operatorname{cone}_{2}=\operatorname{cone}(g_{2},g_{1,2},g_{2,3},g_{1,2,3}), and \operatorname{cone}_{3}=\operatorname{cone}(g_{3},g_{1,3},g_{2,3},g_{1,2,3}). This is illustrated on the right side of Figure [21](https://arxiv.org/html/2603.03425v3#A1.F21 "Figure 21 ‣ Example A.2. ‣ A.2.3 Towards a general construction: 𝒳 in 𝑌 fans in dimensions 2 and 3 ‣ A.2 The 𝒰-polytope for one-loop-integrated amplitudes ‣ Appendix A Appendix: 𝒳 in Y polytopes ‣ Combinatorics of the cosmohedron"). Already in this example, we see that the fan is not simplicial, and therefore the \mathcal{U}-polytope will not be simple.

#### A.2.4 The normal fan of the \mathcal{U}-polytope

Now we are ready to describe the general construction of the \mathcal{U}-fan.

Fan of the one-loop associahedron. The fan of the one-loop associahedron [[AH](https://arxiv.org/html/2603.03425v3#bib.bibx3)] has a ray for each curve of \hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n}, and a face \operatorname{cone}(\mathcal{C}) for each collection of compatible curves, \mathcal{C}\subset\operatorname{Curves}(\hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n}), or equivalently, each subdivision of \hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n}. This is a simplicial fan whose cones are ordered by reverse inclusion:

\operatorname{cone}(\mathcal{C})\subseteq\operatorname{cone}(\mathcal{C}^{\prime})\quad\text{ if and only if }\quad\mathcal{C}\subseteq\mathcal{C}^{\prime}.

Fan of the \mathcal{U} polytope. Now consider a facet \operatorname{cone}(\mathcal{C}) of the one-loop associahedron, corresponding to a triangulation \mathcal{C} of \hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n}. In the polytope, we wish to replace the vertex v_{\mathcal{C}} with a simplex \Delta_{\operatorname{Spokes}(\mathcal{C})}, since the spokes of \mathcal{C} correspond to the edges in the unique loop of the cubic diagram dual to \mathcal{C}. Therefore we need to subdivide the simplicial face \operatorname{cone}(\mathcal{C}) into |\operatorname{Spokes}(\mathcal{C})| maximal cones. If \operatorname{cone}(\mathcal{C}) is regular then the s-th cone \operatorname{cone}(\mathcal{C},s) consists of the points in \operatorname{cone}(\mathcal{C}) whose closest spoke-ray is g_{s}. If it is not, we linearly transform it into a regular cone, and use that subdivision. Explicitly, we obtain new rays

g_{\mathcal{S}}=\sum_{s\in\mathcal{S}}g_{s},\qquad\tilde{g}_{\mathcal{S}}=\sum_{s\in{\mathcal{S}}}\tilde{g}_{s}

for each subset \mathcal{S}\subseteq[n] and more generally, a new cone

\operatorname{cone}(\mathcal{C},\mathcal{S})=\operatorname{cone}(\{g_{ij}\,:\,x_{ij}\in\operatorname{Chords}(\mathcal{C})\}\cup\{g_{T}\,:\,S\subseteq T\subseteq\operatorname{Spokes}(\mathcal{C})\})

for each non-crossing set of curves \mathcal{C} and each non-empty subset \emptyset\neq\mathcal{S}\subseteq\operatorname{Spokes}(\mathcal{C}) – where we allow the exception \mathcal{S}=\emptyset for the sets \mathcal{C} with no spokes.

Note that \operatorname{cone}(\mathcal{C},\mathcal{S}) is a |\operatorname{Chords}(\mathcal{C})|-fold cone over a (|\operatorname{Spokes}(\mathcal{C})|-|S|)–cube. We have

\operatorname{cone}(\mathcal{C},\mathcal{S})\subseteq\operatorname{cone}(\mathcal{C}^{\prime},\mathcal{S}^{\prime})\quad\text{ if and only if }\quad\mathcal{C}\subseteq\mathcal{C}^{\prime}\text{ and }\mathcal{S}\supseteq\mathcal{S}^{\prime},

The blue and red fans glue along the “equatorial complex” of (necessarily lower-dimensional) faces \operatorname{cone}(\mathcal{C},\emptyset) such that \mathcal{C} has no spokes.

#### A.2.5 The \mathcal{U}-polytope

We now construct a geometric realization of the \mathcal{U}-polytope, given by the appropriate chiseling procedure that creates, at each vertex of the underlying Y-polytope, the respective \mathcal{X}-polytope. We understood above how to subdivide the normal fan of the loop associahedron to get the normal fan of the \mathcal{U}-polytope – and therefore we have its new facet normals g_{\mathcal{S}} and \tilde{g}_{\mathcal{S}} for \mathcal{S}\subseteq[n] . However, as we saw for the cosmohedron, it is a subtle matter to determine how deeply we need to chisel in each direction to obtain the correct polytope. This is controlled by the constants \epsilon_{\mathcal{S}} and \tilde{\epsilon}_{\mathcal{S}}.appearing in the new facet inequalities g_{\mathcal{S}}\cdot x\geq\epsilon_{\mathcal{S}} and \tilde{g}_{\mathcal{S}}\cdot x\geq\tilde{\epsilon}_{\mathcal{S}}.

Given a vector \mathsf{c} of positive parameters c_{rs}, c_{r} and \tilde{c}_{r} for 1\leq r,s\leq n, the _one-loop associahedron_ L_{1}(\mathsf{c}) defined in [[AH](https://arxiv.org/html/2603.03425v3#bib.bibx3), [BR25](https://arxiv.org/html/2603.03425v3#bib.bibx13)] is the n-dimensional polytope

\displaystyle L_{1}(\mathsf{c})=\Bigg\{(x,y,\tilde{y})\displaystyle\in\mathbb{R}^{\operatorname{Curves}(\hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n+2})}\,:\,\left[x_{r,s}\geq 0,\,y_{r}\geq 0,\,\tilde{y}_{r}\geq 0\,\right],(13)
\displaystyle\left[\begin{aligned} &x_{r,r}+x_{r+1,r+1}-x_{r+1,r}-y_{r}-\tilde{y}_{r-1}=c_{r+1,r},\\
&x_{r,s+1}+x_{r-1,s}-x_{r,s}-x_{r-1,s-1}=c_{rs}\text{ for }r\neq\{s+1,s\},\\
&y_{r}+\tilde{y}_{r}-x_{r+1,r+1}=c_{r+1}.\end{aligned}\right]\Bigg\}

where we set to 0 any variables that are undefined because their subscripts are outside of the prescribed ranges. In this equation, we are writing vectors in \mathbb{R}^{\operatorname{Curves}(\hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n+2})} in the form (x,y,\tilde{y})\in\mathbb{R}^{\operatorname{Chords}(\hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n+2})}\times\mathbb{R}^{\operatorname{Spokes}(\hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n+2})}\times\mathbb{R}^{\widetilde{\operatorname{Spokes}}(\hbox to11.22pt{\vbox to10.69pt{\pgfpicture\makeatletter\hbox{\hskip 5.61198pt\lower-4.8037pt\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} \lxSVG@begingroup@{stroke=#000000} \lxSVG@begingroup@{fill=#000000} \lxSVG@setlinewidth{\the\pgflinewidth}\lxSVG@begingroup@{stroke-width=0.4pt} \lx@inpgf@ignorespaces\nullfont\hbox to0.0pt{\lxSVG@begingroup@{_scopebegin=1} {{\lx@inpgf@ignorespaces}}{
{}{{}}{}
{}{}
{}{}
{}{}
{}{}
{\lx@inpgf@ignorespaces}{}\lxSVG@stroke\lxSVG@drawpath@unclipped{M 0 7.87 L -7.49 2.43 L -4.63 -6.37 L 4.63 -6.37 L 7.49 2.43 Z}{fill:none} \lx@inpgf@ignorespaces
{}{{}}{}{{{}}{\lx@inpgf@ignorespaces}{}{}{}{}{}{}{}}\lxSVG@fill\lxSVG@drawpath@unclipped{M 0 0 M 1.18 0 C 1.18 0.65 0.65 1.18 0 1.18 C -0.65 1.18 -1.18 0.65 -1.18 0 C -1.18 -0.65 -0.65 -1.18 0 -1.18 C 0.65 -1.18 1.18 -0.65 1.18 0 Z M 0 0}{stroke:none} \lx@inpgf@ignorespaces
}
\lxSVG@closescope {\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}{\lx@inpgf@ignorespaces}\hss}\lxSVG@discardpath\lxSVG@closescope \hss}}\lxSVG@closescope\endpgfpicture}}_{n+2})} where x=(x_{r,s})_{1\leq r,s\leq n}, y=(y_{r})_{1\leq r\leq n}, and \tilde{y}=(\tilde{y}_{r})_{1\leq r\leq n}.

To construct our \mathcal{U}-polytope, we need to introduce new facet inequalities of the form:

\sum_{k\in\mathcal{S}}y_{k}\geq\epsilon_{\mathcal{S}},\quad\sum_{k\in\mathcal{S}}\tilde{y}_{k}\geq\tilde{\epsilon}_{\mathcal{S}}\qquad\text{ for }\emptyset\neq\mathcal{S}\subseteq[n](14)

for the possible subsets of spokes. A first condition on the \epsilon s is that they should be small enough in relation to the c s that the chiseling will not remove any faces of L_{1}(\mathsf{c}). The parameters must also satisfy the equalities

\epsilon_{\mathcal{S}_{1}}+\epsilon_{\mathcal{S}_{2}}=\epsilon_{\mathcal{S}_{1}\cup\mathcal{S}_{2}}+\epsilon_{\mathcal{S}_{1}\cap\mathcal{S}_{2}}\text{ for any subsets }\mathcal{S}_{1},\mathcal{S}_{2}\text{ with }\mathcal{S}_{1}\cap\mathcal{S}_{2}\neq\emptyset(15)

since the corresponding rays of the normal fan satisfy g_{\mathcal{S}_{1}}+g_{\mathcal{S}_{2}}=g_{\mathcal{S}_{1}\cup\mathcal{S}_{2}}+g_{\mathcal{S}_{1}\cap\mathcal{S}_{2}}, and they lie in the same cone \operatorname{cone}([n],\{y_{s}\}) of the fan for any s\in\mathcal{S}_{1}\cap\mathcal{S}_{2}. In addition, the \epsilon s must satisfy the _wall-crossing inequalities_

\epsilon_{i,j}>\epsilon_{i}+\epsilon_{j}\qquad\text{ for }i\neq j(16)

since the corresponding rays satisfy g_{i,j}=g_{i}+g_{j}, and \operatorname{cone}(g_{i},g_{i,j}) and \operatorname{cone}(g_{j},g_{i,j}) are distinct cones in the fan. The analogous conditions hold for the \tilde{\epsilon}’s. One may verify that these are all the necessary and sufficient conditions to obtain the desired \mathcal{U}-polytope.

The solutions to the system of equations ([15](https://arxiv.org/html/2603.03425v3#A1.E15 "In A.2.5 The 𝒰-polytope ‣ A.2 The 𝒰-polytope for one-loop-integrated amplitudes ‣ Appendix A Appendix: 𝒳 in Y polytopes ‣ Combinatorics of the cosmohedron")) are of the form

\epsilon_{\{s_{1},\ldots,s_{k}\}}=\epsilon_{s_{1}}+\cdots+\epsilon_{s_{k}}+(k-1)\epsilon(17)

for any choice of \epsilon_{1},\ldots,\epsilon_{n},\epsilon, and the inequalities ([16](https://arxiv.org/html/2603.03425v3#A1.E16 "In A.2.5 The 𝒰-polytope ‣ A.2 The 𝒰-polytope for one-loop-integrated amplitudes ‣ Appendix A Appendix: 𝒳 in Y polytopes ‣ Combinatorics of the cosmohedron")) are satisfied when \epsilon>0. The same analysis holds for the \tilde{\epsilon}’s. Thus a particularly simple choice of parameters that cuts out a \mathcal{U}-polytope is

\epsilon_{S}=\tilde{\epsilon}_{S}=(|S|-1)\epsilon\,.(18)

for small \epsilon>0. Summarizing, the \mathcal{U}-polytope is

\mathcal{U}(\mathsf{c})=\{(x,y,\tilde{y})\in L_{1}(\mathsf{c})\,:\,\sum_{k\in\mathcal{S}}y_{k}\geq(|S|-1)\epsilon,\quad\sum_{k\in\mathcal{S}}\tilde{y}_{k}\geq(|S|-1)\epsilon\text{ for }\emptyset\neq\mathcal{S}\subseteq[n]\}

Figure 22: (a) Y: Two-dimensional one-loop associahedron. (b) \mathcal{X} in Y: Two-dimensional one-loop \mathcal{U}-polytope. 

###### Example A.3.

For n=2, the one-loop associahedron is the intersection of the subspace

\displaystyle y_{1}+\tilde{y}_{1}-x_{2,2}=c_{2}\,,\quad y_{2}+\tilde{y}_{2}-x_{1,1}=c_{1}\,,
\displaystyle x_{1,1}+x_{2,2}-y_{2}-\tilde{y}_{1}=c_{1,2}\,,\quad x_{1,1}+x_{2,2}-y_{1}-\tilde{y}_{2}=c_{2,1}\,,(19)

with the positive orthant x_{1,1}\geq 0,\,x_{2,2}\geq 0,\,y_{1}\geq 0,\,y_{2}\geq 0. To obtain the \mathcal{U}-polytope, we impose the stronger inequalities

\displaystyle x_{1,1}\geq 0\,,\,\,x_{2,2}\geq 0\,,\,\,y_{1}\geq\epsilon_{1}\,,\,\,y_{2}\geq\epsilon_{2}\,,\,\,y_{1}+y_{2}\geq\epsilon_{1}+\epsilon_{2}+\epsilon\,.(20)

for any \epsilon_{1},\epsilon_{2},\epsilon>0 and the analogous inequalities for \tilde{y} for any \tilde{\epsilon}_{1},\tilde{\epsilon}_{2},\tilde{\epsilon} with \epsilon,\tilde{\epsilon}>0. Using the parameterization \epsilon_{1}=\epsilon_{2}=0 (and similarly for \tilde{\epsilon}) and \epsilon=\tilde{\epsilon} we obtain the correct embedding as shown on the right of Figure [22](https://arxiv.org/html/2603.03425v3#A1.F22 "Figure 22 ‣ A.2.5 The 𝒰-polytope ‣ A.2 The 𝒰-polytope for one-loop-integrated amplitudes ‣ Appendix A Appendix: 𝒳 in Y polytopes ‣ Combinatorics of the cosmohedron").

Figure 23: (a) Y: Three-dimensional one-loop associahedron. (b) \mathcal{X} in Y: Three dimensional one-loop \mathcal{U}-polytope. 

Similarly, for n=3, we intersect the subspace

\displaystyle y_{1}+\tilde{y}_{1}-x_{2,2}=c_{2}\,,\quad y_{2}+\tilde{y}_{2}-x_{3,3}=c_{3}\,,\quad y_{3}+\tilde{y}_{3}-x_{1,1}=c_{1}\,,
\displaystyle x_{1,1}+x_{2,2}+x_{2,1}-y_{1}-\tilde{y}_{3}=c_{2,1}\,,\,\,\,x_{2,2}+x_{3,3}+x_{3,2}-y_{2}-\tilde{y}_{1}=c_{3,2}
\displaystyle x_{3,3}+x_{1,1}+x_{1,3}-y_{3}-\tilde{y}_{2}=c_{1,3}\,,\quad x_{1,3}+x_{2,1}-x_{1,1}=c_{2,3}\,,
\displaystyle x_{2,1}+x_{3,2}-x_{2,2}=c_{3,1}\,,\quad x_{3,2}+x_{1,3}-x_{3,3}=c_{1,2}\,,(21)

with the inequalities

\displaystyle x_{i,i}\geq 0\,,\quad x_{2,1}\geq 0\,,\quad x_{3,2}\geq 0\,,\quad x_{1,3}\geq 0,
\displaystyle y_{i}\geq\epsilon_{i}\,,\quad y_{i}+y_{j}\geq\epsilon_{i}+\epsilon_{j}+\epsilon\,,\quad y_{1}+y_{2}+y_{3}\geq\epsilon_{1}+\epsilon_{2}+\epsilon_{3}+2\epsilon(22)

for \epsilon_{1},\epsilon_{2},\epsilon_{3},\epsilon>0, and analogous inequalities for the \tilde{y}s.

With \epsilon_{i}=\tilde{\epsilon}_{i}=0 and \epsilon=\tilde{\epsilon}=1, we produce the correct embedding seen in Figure [23](https://arxiv.org/html/2603.03425v3#A1.F23 "Figure 23 ‣ Example A.3. ‣ A.2.5 The 𝒰-polytope ‣ A.2 The 𝒰-polytope for one-loop-integrated amplitudes ‣ Appendix A Appendix: 𝒳 in Y polytopes ‣ Combinatorics of the cosmohedron").

## References

*   [AA23] Marcelo Aguiar and Federico Ardila. Hopf monoids and generalized permutahedra, volume 289. American Mathematical Society, 2023. 
*   [ACEP20] Federico Ardila, Federico Castillo, Christopher Eur, and Alexander Postnikov. Coxeter submodular functions and deformations of Coxeter permutahedra. Advances in Mathematics, 365:107039, 2020. 
*   [AH] Nima Arkani-Hamed. All-Loop Scattering as a Counting Problem. Amplitudes 2022 Talk. 
*   [AHBHY18a] Nima Arkani-Hamed, Yuntao Bai, Song He, and Gongwang Yan. Scattering Forms and the Positive Geometry of Kinematics, Color and the Worldsheet. JHEP, 05:096, 2018. 
*   [AHBHY18b] Nima Arkani-Hamed, Yuntao Bai, Song He, and Gongwang Yan. Scattering forms and the positive geometry of kinematics, color and the worldsheet. Journal of High Energy Physics, 2018(5):1–78, 2018. 
*   [AHCD+24] Nima Arkani-Hamed, Qu Cao, Jin Dong, Carolina Figueiredo, and Song He. Hidden zeros for particle/string amplitudes and the unity of colored scalars, pions and gluons. JHEP, 10:231, 2024. 
*   [AHFS+25] N.Arkani-Hamed, H.Frost, G.Salvatori, P-G. Plamondon, and H.Thomas. All loop scattering as a counting problem. JHEP, 08:194, 2025. 
*   [AHFV25] Nima Arkani-Hamed, Carolina Figueiredo, and Francisco Vazão. Cosmohedra. JHEP, 11:029, 2025. 
*   [AHHM22] Nima Arkani-Hamed, Aaron Hillman, and Sebastian Mizera. Feynman polytopes and the tropical geometry of UV and IR divergences. Phys. Rev. D, 105(12):125013, 2022. 
*   [AHHST22] Nima Arkani-Hamed, Song He, Giulio Salvatori, and Hugh Thomas. Causal diamonds, cluster polytopes and scattering amplitudes. JHEP, 11:049, 2022. 
*   [ARW06] Federico Ardila, Victor Reiner, and Lauren Williams. Bergman complexes, coxeter arrangements, and graph associahedra. Séminaire Lotharingien de Combinatoire, 54:B54Aj, 2006. 
*   [BCL+22] Luciana Basualdo Bonatto, Safia Chettih, Abigail Linton, Sophie Raynor, Marcy Robertson, and Nathalie Wahl. An infinity operad of normalized cacti. Topology and its Applications, 316:108107, 2022. 
*   [BR25] Jeffrey V. Backus and Laurentiu Rodina. Emergence of Unitarity and Locality from Hidden Zeros at One-Loop Order. Phys. Rev. Lett., 135(13):131601, 2025. 
*   [BS13] Alicia Boole Stott. Geometrical deduction of semiregular from regular polytopes and space fillings, volume 11. J. Müller, 1913. 
*   [CD06] Michael Carr and Satyan L Devadoss. Coxeter complexes and graph-associahedra. Topology and its Applications, 153(12):2155–2168, 2006. 
*   [CL22] Federico Castillo and Fu Liu. Deformation cones of nested braid fans. International Mathematics Research Notices, 2022(3):1973–2026, 2022. 
*   [CL23] Federico Castillo and Fu Liu. The permuto-associahedron revisited. European Journal of Combinatorics, 110:103706, 2023. 
*   [CLS24] David A Cox, John B Little, and Henry K Schenck. Toric varieties, volume 124. American Mathematical Society, 2024. 
*   [CSZ15] Cesar Ceballos, Francisco Santos, and Günter M Ziegler. Many non-equivalent realizations of the associahedron. Combinatorica, 35(5):513–551, 2015. 
*   [Dev09] Satyan L Devadoss. A realization of graph associahedra. Discrete Mathematics, 309(1):271–276, 2009. 
*   [Ear24] Nick Early. Generalized permutohedra in the kinematic space. Journal of High Energy Physics, 2024(6):1–24, 2024. 
*   [Edm70] Jack Edmonds. Submodular functions, matroids, and certain polyhedra. In Combinatorial Structures and their Applications (Proc. Calgary Internat. Conf., Calgary, Alta., 1969), pages 69–87. Springer, 1970. 
*   [FV26] Carolina Figueiredo and Francisco Vazão. Correlator polytopes. Phys. Rev. D, 113(2):025005, 2026. 
*   [Gai15] Giovanni Gaiffi. Permutonestohedra. Journal of Algebraic Combinatorics, 41(1):125–155, 2015. 
*   [Kap93] Mikhail M Kapranov. The permutoassociahedron, mac lane’s coherence theorem and asymptotic zones for the kz equation. Journal of Pure and Applied Algebra, 85(2):119–142, 1993. 
*   [KI14] Vaclav Kotesovic and OEIS Foundation Inc. Entry E177384, The on-line encyclopedia of integer sequences, 2014. 
*   [Knu97] Donald E Knuth. The Art of Computer Programming: Fundamental Algorithms, Volume 1. Addison-Wesley Professional, 1997. 
*   [LA22] Guillaume Laplante-Anfossi. The diagonal of the operahedra. Advances in Mathematics, 405:108494, 2022. 
*   [Lod04] Jean-Louis Loday. Realization of the stasheff polytope. Archiv der Mathematik, 83(3):267–278, 2004. 
*   [Mel21] Stephen Melczer. An Invitation to Analytic Combinatorics. Springer, 2021. 
*   [MSS02] Martin Markl, Steve Shnider, and Jim Stasheff. Operads in algebra, topology and physics. Mathematical surveys and monographs, 96, 2002. 
*   [Pos09] Alexander Postnikov. Permutohedra, associahedra, and beyond. International Mathematics Research Notices, 2009(6):1026–1106, 2009. 
*   [PPP23] Arnau Padrol, Vincent Pilaud, and Germain Poullot. Deformation cones of graph associahedra and nestohedra. European Journal of Combinatorics, 107:103594, 2023. 
*   [Sch13] Pieter Hendrik Schoute. Analytical treatment of the polytopes regularly derived from the regular polytopes, volume 11. J. Müller, 1913. 
*   [Sta63] James Dillon Stasheff. Homotopy associativity of h-spaces. ii. Transactions of the American Mathematical Society, 108(2):293–312, 1963. 
*   [Sta97] Richard P Stanley. Hipparchus, plutarch, schröder, and hough. The American Mathematical Monthly, 104(4):344–350, 1997. 
*   [Sta99] Richard P Stanley. Enumerative combinatorics volume 2. Cambridge studies in advanced mathematics, 1999. 
*   [Tam54] Dov Tamari. Monoïdes préordonnés et chaînes de malcev. Bulletin de la Société mathématique de France, 82:53–96, 1954.
