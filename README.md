# Improving Single and Multi-View Blockmodelling by Algebraic Simplification

This contains the work done as a part of my 2019 summer internship at Monash and University of Melbourne under the guidance of [Prof. Peter Stuckey](https://scholar.google.com/citations?user=tvFekxwAAAAJ&hl=en) and [Prof. Jefrey Chan](https://scholar.google.com.au/citations?user=I9mDeiYAAAAJ&hl=en)

Abstract:
Blockmodelling is an important technique in social network analysis for discovering the latent structures and groupings in graphs. State-of-the-art approaches approximate the graph using matrix factorisation, which can discover both the latent graph structures and vertex groupings. However, factorisation is a one-way approximation, in that it only approximates the graph with a lossy model that removes the background noise. Traditional Blockmodelling methods rely on an alternating 2-step optimization that involves iteratively updating the matrix representing membership while fixing the matrix representing the graph's underlying structure, and then updating the structure matrix while keeping the membership matrix fixed. We propose a single step optimization method, which uses algebraic simplifi-cation to directly update the lower dimensional, latent structure representation. This helps improve both the convergence and accuracy of blockmodelling. We also show that this approach can solve multi-view blockmodelling problems, involving multiple graphs over the same vertices. We use real datasets to show that our approach has much higher accuracy and comparable running times to competing approaches.

Published in **2020 International Joint Conference on Neural Networks (IJCNN), Glasgow(UK)**


 ([Link to paper](https://ieeexplore.ieee.org/document/9207065))

