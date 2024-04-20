# HAPAGM
Heurística de Particionamento de Árvore Geradora Mínima
Esta heurística corresponde corresponde a uma variante de outro algoritmo proposto por (Brito e Montenegro, 2010) e tem,
por base, a construção e particionamento de uma árvore geradora mínima (AGM) T. A partir de T são gerados todas
os agrupamentos espaciais (subárvores),  considerando todas as possíveis remoções de k-1 arestas dentre n-1 arestas disponíveis em T. 
O total de possíveis agrupamentos espaciais obtidos a partir da remoção de (k-1) arestas é dado por (n-1)! / [(n-1-k+1)! (k-1)!]
Para cada agrupamento espacial produzido avalia-se uma restrição de capacidade (mínimo de objetos por grupo) e uma função objetivo,
