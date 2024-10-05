PAUTET Florence

# PHYG TP: comparasion de séquences sans alignements

## 1: Notions générales

L'objectif de ce TP est de comparer des séquences sans avoir à les aligner. Pour cela, nous utilisons les k-mers (mots de longueur $k$) qu'elles contiennent. Chaque séquence est alors décrite par le multiset de ses k-mers.

Pour comparer deux séquences de multisets $A$ et $B$, nous utilisons le coefficient de Jaccard:

$$J(A,B) = \frac{|A\cap B|}{|A\cup B|}$$

## 2: Présentation du code

### Encodage des lettres
Pour les séquences ADN, l'alphabet utilisé contient 4 lettres: {A,C,G,T}. Elles peuvent donc être codées sur 2 bits. Une façon efficace de les encoder est d'utiliser les 2ème et 3ème bit de leur code ASCII.

Pour passer du code $c$ d'un nucléotide à celui $\bar{c}$ de son complémentaire, on peut effectuer un _ou exclusif_ avec un masque adapté ${(0B10)}$.

| lettre | code ASCII | code | code du complémentaire |
|:------:|:----------:|:----:|:----------------------:|
| A | 01000**00**1 | 00 | 10 |
| C | 01000**01**1 | 01 | 11 |
| T | 01010**10**0 | 10 | 00 |
| G | 01000**11**1 | 11 | 01 |


Ainsi un kmer $a_1, \ldots, a_k$ est encodé sur $2k$ bits par la séquence $c_1,\ldots, c_k$.  
De plus, sa séquence inverse complémentaire est codé par $\bar{c_k},\ldots, \bar{c_1}$.

### Énumération des k-mers

Pour énumérer les k-mers d'une séquence, nous l'itérons et construisons successivement chaque k-mer et son complémentaire inverse. Pour passer de $c_{i},\ldots, c_{i+k-1}$ à  $c_{i+1},\ldots, c_{i+k}$, il suffit de décaler de deux bits vers la gauche puis d'ajouter $c_{i+k}$ (à droite). De même pour le complémentaire inverse, il faut décaler de deux bits vers la droite puis ajouter ${\bar{c}_{i+k}<<2(k-1)}$ (à gauche).

Nous ne faisons pas la différence entre une séquence et son inverse complémentaire. Pour chaque kmer, nous retournons donc uniquement la plus petite valeur entre la sienne et celle de son inverse complémentaire.

### Construction du dictionnaire

Chaque organisme est représentée par un dictionaire {kmer: nombre d'occurences}. Ce dictionnaire est construit par itérations. Dans notre cas, les organismes sont décrits par des listes de plusieurs séquences. Le multi-set contient alors le nombre d'occurences total.

### Calcul du coefficient de Jaccard

Nous précalculons les dictionnaires représentatifs de chaque organisme. Ensuite, pour chaque paire d'organismes, nous calculons la distance de Jaccard à l'aide des dictionnaires.

#### Coefficient de Jaccard:
$$ \begin{align*}
J(A,B) & =  \frac{|A\cap B|}{|A\cup B|} \\
       & =  \frac{|A\cap B|}{|A| + |B| - |A\cap B|}
\end{align*}$$

$|A\cap B|$ est calculé en itérant un dictionnaire et en vérifiant si la clé est présente dans l'autre dictionnaire. Si oui, on prend le minimum des deux valeurs.

$|A|$ et $|B|$ peuvent être précalculés au moment de la création des dictionaire.

### Bonus: Tests unitaires
Toutes les fonctions d'encodage ont été vérifiées via des test untaires présents dans `test.py`.

## 3: Analyse des résultas

$$
\begin{array}{l|ccccc}
    & \text{GCA\\_000005845.2} & \text{GCA\\_000013265.1} & \text{GCA\\_030271835.1} & \text{GCA\\_000069965.1} & \text{GCA\\_000008865.2} \\
    \hline
    \text{GCA\\_000005845.2} & 1 & 0.3410 & 0.0026 & 0.0026 & 0.4365 \\
    \text{GCA\\_000013265.1} & 0.3410 & 1 & 0.0024 & 0.0024 & 0.3071 \\
    \text{GCA\\_030271835.1} & 0.0026 & 0.0024 & 1 & 0.0311 & 0.0023 \\
    \text{GCA\\_000069965.1} & 0.0026 & 0.0024 & 0.0311 & 1 & 0.0023 \\
    \text{GCA\\_000008865.2} & 0.4365 & 0.3071 & 0.0023 & 0.0023 & 1 \\
\end{array}
$$

Ceci est une matrice de similarité. Plus la valeur est proche de 1, plus les organismes sont proches.

En appliquant un algorithme de type UPGMA, on obtient un groupe de 3 organismes relativement similaires (GCA_000005845.2, GCA_000013265.1 et GCA_000008865.2) et un groupe de deux un peu plus éloignés. 

Cela fait sens, puisque qu'en regardant les organismes correspondant, on constate que les séquences les plus proches appartiennent à _Escherichia coli_, tandis que les deux autres organismes sont deux bactéries du genre _Proteus_: _Proteus appendicitidis_ et _Proteus mirabilis_.

