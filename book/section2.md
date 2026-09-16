# 2. Théorème Spectral

Le Théorème Spectral permet de représenter l'application linéaire $T : \mathbb{R}^n \to \mathbb{R}^n$ associée à $A \in \mathbb{R}^{n \times n}$ une matrice symétrique réelle si et seulement si $A$ peut être exprimée sous la forme :

$$
A = QDQ^{\top}
$$

où $Q$ est une matrice orthogonale dont les colonnes $\{u_1, \dots, u_n\}$ sont les vecteurs propres de $A$. Ces colonnes forment une base orthonormale pour $\mathbb{R}^n$. Puisque $Q$ est orthogonal, son inverse est égal à sa transposée, ${Q^{-1} = Q^{\top}}$.

$D$ est une matrice diagonale dont les éléments diagonaux sont les valeurs propres de $A$, correspondant aux vecteurs propres dans $Q$ dans le même ordre.

Autrement dit, on peut écrire la matrice $A$ simplement comme une somme de ses valeurs propres et vecteurs propres :

$$
A = \lambda_1 u_1 u_1^\top + \lambda_2 u_2 u_2^\top + \dots + \lambda_n u_n u_n^\top = \sum_{i=1}^{n} \lambda_i u_i u_i^\top.
$$

```{admonition} Exemple
:class: seealso
Considérons la matrice symétrique $A\in \mathbb{R}^{3 \times 3}$ suivante :

$$
A =
\begin{pmatrix}
    2 & 0 & 0 \\
    0 & 3 & 1 \\
    0 & 1 & 3
\end{pmatrix}.
$$
```

Tout d'abord on détermine $\operatorname{spec}(A)$.
Le polynôme caractéristique est

$$
\det(A - \lambda I) =
\det
\begin{pmatrix}
    2-\lambda & 0         & 0         \\
    0         & 3-\lambda & 1         \\
    0         & 1         & 3-\lambda
\end{pmatrix} = 0.
$$

On obtient :

$$
(2-\lambda)((3-\lambda)^2 - 1) = 0 \implies (2-\lambda)(2-\lambda)(4-\lambda) = 0 \implies (2-\lambda)^2 (4-\lambda) = 0.
$$

Ainsi, l'ensemble des valeurs propres est $\operatorname{spec}(A)$ = {2, 4}.
Ensuite, on trouve les vecteurs $v$ non-nuls qui satisfait $(A - \lambda I)v = 0$.

Pour $\lambda = 2$:

$$
(A - 2I)v = 0 \implies
\begin{pmatrix} 0 & 0 & 0 \\ 0 & 1 & 1 \\ 0 & 1 & 1 \end{pmatrix}
\begin{pmatrix}x\\y\\z\end{pmatrix} =
\begin{pmatrix}0\\0\\0\end{pmatrix} \implies y + z = 0.
$$

Les vecteurs de base sont : $v_1 = \begin{pmatrix}1\\0\\0\end{pmatrix}, \quad v_2 = \begin{pmatrix}0\\1\\-1\end{pmatrix}.$

Pour $\lambda = 4$:

$$
(A - 4I)v = 0 \implies
\begin{pmatrix} -2 & 0 & 0 \\ 0 & -1 & 1 \\ 0 & 1 & -1 \end{pmatrix}
\begin{pmatrix}x\\y\\z\end{pmatrix} =
\begin{pmatrix}0\\0\\0\end{pmatrix} \implies y=z.
$$

Le vecteur de base est : $v_3 = \begin{pmatrix}0\\1\\1\end{pmatrix}.$

On va écrire une matrice orthogonale $Q$. Cela signifie que les vecteurs propres que nous plaçons dans $Q$ peuvent être et orthogonaux entre eux. Ainsi, si un vecteur propre n'est pas de longueur 1, nous le normalisons en le divisant par sa norme $\|v\|$ : Nous normalisons les vecteurs $v$ en $u = v/\|v\|$:

$$
u_1 = \frac{v_1}{\|v_1\|} = \frac{1}{1} \begin{pmatrix}1\\0\\0\end{pmatrix} = \begin{pmatrix}1\\0\\0\end{pmatrix}
$$

$$
u_2 = \frac{v_2}{\|v_2\|} = \frac{1}{\sqrt{0^2+1^2+(-1)^2}} \begin{pmatrix}0\\1\\-1\end{pmatrix} = \begin{pmatrix}0\\1/\sqrt{2}\\-1/\sqrt{2}\end{pmatrix}
$$

$$
u_3 = \frac{v_3}{\|v_3\|} = \frac{1}{\sqrt{0^2+1^2+1^2}} \begin{pmatrix}0\\1\\1\end{pmatrix} = \begin{pmatrix}0\\1/\sqrt{2}\\1/\sqrt{2}\end{pmatrix}
$$

$Q$ est la matrice orthogonale formée des vecteurs propres normalisés $\{u_1, u_2, u_3\}$:

$$
Q = [u_1 \; u_2 \; u_3] =
\begin{pmatrix}
    1 & 0           & 0          \\
    0 & 1/\sqrt{2}  & 1/\sqrt{2} \\
    0 & -1/\sqrt{2} & 1/\sqrt{2}
\end{pmatrix}.
$$

$D$ est la matrice diagonale des valeurs propres correspondantes:

$$
D = \text{diag}(2, 2, 4) =
\begin{pmatrix}
    2 & 0 & 0 \\
    0 & 2 & 0 \\
    0 & 0 & 4
\end{pmatrix}.
$$

Le Théorème Spectral est vérifié :

$$
\boxed{A = Q D Q^{\top}}
$$

```{admonition} Lemme (Valeurs propres réelles)
:class: note
Toutes les valeurs propres d'une matrice symétrique réelle sont réelles.
```

```{admonition} Démonstration
Soit $\lambda \in \mathbb{C}$ une valeur propre de la matrice réelle symétrique $A$ et $x \in \mathbb{C}^n$ un vecteur propre associé non nul ($x \neq 0$), tel que $Ax = \lambda x$.

On considère $\mathbb{C}^n$ muni du produit scalaire hermitien usuel noté $\langle \cdot, \cdot \rangle$.
Puisque $A$ est réelle et symétrique ($A = A^\top$), on a $A = \overline{A}^\top = A^*$.

Par la propriété de l'adjoint dans un produit scalaire ($\langle Ax, y \rangle = \langle x, A^* y \rangle$), nous obtenons l'égalité suivante :

$$
\langle Ax, x \rangle = \langle x, A^* x \rangle = \langle x, Ax \rangle.
$$

En remplaçant $A x = \lambda x$, on obtient :

$$
\langle \lambda x, x \rangle = \langle x, \lambda x \rangle
$$

En utilisant les propriétés de linéarité du produit scalaire :

$$
\lambda \langle x, x \rangle = \bar{\lambda} \langle x, x \rangle
$$

Comme $x \neq 0$, on a $\langle x, x \rangle = \|x\|^2 > 0$. Nous pouvons donc simplifier l'équation :

$$
\lambda = \bar{\lambda}
$$

Ceci implique que la partie imaginaire de $\lambda$ est nulle, donc $\lambda \in \mathbb{R}$.
```

```{admonition} Lemme (Vecteurs propres orthogonaux)
:class: note
Les vecteurs propres correspondant à des valeurs propres distinctes d'une matrice symétrique sont orthogonaux.
```

```{admonition} Démonstration
Soient $\lambda_1$ et $\lambda_2$ deux valeurs propres distinctes correspondant respectivement aux vecteurs propres $x_1$ et $x_2$, tels que $Ax_1 = \lambda_1 x_1$ et $Ax_2 = \lambda_2 x_2$.

Nous allons établir deux égalités pour comparer les produits scalaires.
Prenons la première équation $Ax_1 = \lambda_1 x_1$. Appliquons d'abord la transposée des deux côtés, puis multiplions à droite par le vecteur $x_2$ :

$$
(A x_1)^\top = (\lambda_1 x_1)^\top \implies x_1^\top A^\top = \lambda_1 x_1^\top
$$

Puisque $A$ est symétrique ($A = A^\top$), on remplace $A^\top$ par $A$. En multipliant ensuite à droite par $x_2$, on obtient :

$$
x_1^\top A x_2 = \lambda_1 x_1^\top x_2
$$

Partons maintenant de la deuxième équation $Ax_2 = \lambda_2 x_2$. Cette fois, multiplions directement à gauche par le vecteur transposé $x_1^\top$ :

$$
x_1^\top (A x_2) = x_1^\top (\lambda_2 x_2) \implies x_1^\top A x_2 = \lambda_2 x_1^\top x_2
$$

En égalisant ces deux expressions, nous obtenons :

$$
\lambda_1 x_1^\top x_2 = \lambda_2 x_1^\top x_2
$$

Ce qui équivaut à :

$$
(\lambda_1 - \lambda_2) x_1^\top x_2 = 0
$$

Puisque les valeurs propres sont distinctes ($\lambda_1 \neq \lambda_2$), on a nécessairement $\lambda_1 - \lambda_2 \neq 0$. Par conséquent, le produit scalaire doit être nul :

$$
x_1^\top x_2 = 0
$$

Cela prouve que les vecteurs propres $x_1$ et $x_2$ sont orthogonaux.
```

```{admonition} Lemme (Indépendance des Vecteurs Propres)
:class: note
Si une matrice $A \in \mathbb{R}^{n \times n}$ possède $k$ valeurs propres distinctes, alors tout ensemble de $k$ vecteurs propres correspondants (non nuls) est linéairement indépendant.
```

```{admonition} Démonstration
Nous prouvons ce théorème par récurrence.

**Initialisation :** D'abord, montrons que deux vecteurs propres correspondant à des valeurs propres distinctes sont linéairement indépendants. Soient $v_1$ et $v_2$ des vecteurs propres associés respectivement aux valeurs propres distinctes $\lambda_1$ et $\lambda_2$. Supposons, par l'absurde, que $v_1$ et $v_2$ sont linéairement dépendants. Alors il existe deux scalaires $\alpha_1, \alpha_2$ (non tous les deux nuls) tels que:

$$
\alpha_1 v_1 + \alpha_2 v_2 = 0
$$ (eq:indep1)

En multipliant ({eq}`eq:indep1`) à gauche par $A$, on obtient :

$$
\alpha_1 \lambda_1 v_1 + \alpha_2 \lambda_2 v_2 = 0
$$ (eq:indep2)

De même, en multipliant ({eq}`eq:indep1`) par $\lambda_2$, on obtient :

$$
\alpha_1 \lambda_2 v_1 + \alpha_2 \lambda_2 v_2 = 0
$$ (eq:indep3)

En soustrayant l'équation ({eq}`eq:indep2`) de l'équation ({eq}`eq:indep3`), on trouve :

$$
\alpha_1 (\lambda_2 - \lambda_1) v_1 = 0.
$$

Puisque $\lambda_2 \neq \lambda_1$ et $v_1 \neq 0$, nous devons avoir $ \alpha_1 = 0$. Comme $v_2 \neq 0$, en remplaçant $\alpha_1=0$ dans ({eq}`eq:indep1`), on obtient $\alpha_2=0$, ce qui mène à une contradiction. Ainsi, $v_1$ et $v_2$ sont linéairement indépendants.

**Hérédité :** Supposons maintenant que tout ensemble de $j < k$ vecteurs propres correspondant à des valeurs propres distinctes est linéairement indépendant. Nous voulons montrer que $j+1$ vecteurs propres sont aussi linéairement indépendants.

Soient $v_1, v_2, \dots, v_j$ des vecteurs propres linéairement indépendants correspondant à des valeurs propres distinctes $\lambda_1, \lambda_2, \dots, \lambda_j$. Supposons par l'absurde qu'un vecteur propre supplémentaire $v_{j+1}$, correspondant à une valeur propre différente $\lambda_{j+1}$, est linéairement dépendant de $v_1, \dots, v_j$. Alors il existe des scalaires $\alpha_1, \alpha_2, \dots, \alpha_j$ non tous nuls tels que :

$$
v_{j+1} = \alpha_1 v_1 + \alpha_2 v_2 + \dots + \alpha_j v_j
$$ (eq:indep4)

En multipliant ({eq}`eq:indep4`) à gauche par $A$ :

$$
\lambda_{j+1} v_{j+1} = \alpha_1 \lambda_1 v_1 + \alpha_2 \lambda_2 v_2 + \dots + \alpha_j \lambda_j v_j
$$ (eq:indep5)

En multipliant ({eq}`eq:indep4`) par le scalaire $\lambda_{j+1}$ :

$$
\lambda_{j+1} v_{j+1} = \alpha_1 \lambda_{j+1} v_1 + \alpha_2 \lambda_{j+1} v_2 + \dots + \alpha_j \lambda_{j+1} v_j
$$ (eq:indep6)

En soustrayant ces deux équations, on obtient :

$$
\alpha_1 (\lambda_{j+1} - \lambda_1) v_1 + \alpha_2 (\lambda_{j+1} - \lambda_2) v_2 + \dots + \alpha_j (\lambda_{j+1} - \lambda_j) v_j = 0.
$$

D'après l'hypothèse de récurrence, les vecteurs $v_i$ (pour $i \le j$) sont linéairement indépendants. Par conséquent, tous les coefficients doivent être nuls :

$$
\alpha_i (\lambda_{j+1} - \lambda_i) = 0 \quad \text{pour tout } i \in \{1, \dots, j\}.
$$

Comme les valeurs propres sont distinctes ($\lambda_{j+1} \neq \lambda_i$), on a nécessairement $\alpha_i = 0$ pour tout $i$. Cela contredit l'hypothèse que les vecteurs sont dépendants. Ainsi, les vecteurs propres $v_1, \dots, v_{j+1}$ sont linéairement indépendants.

Par récurrence, tout ensemble de $k$ vecteurs propres correspondant à $k$ valeurs propres distinctes est linéairement indépendant.
```

```{admonition} Exemple
:class: seealso
Considérons la matrice $3 \times 3$ définie par :

$$
A = \begin{pmatrix}
    4 & 0  & 0  \\
    0 & 2  & -2 \\
    0 & -2 & 2
\end{pmatrix}
$$
```

Elle possède trois valeurs propres (calculons $\det(A-\lambda I) = (4-\lambda)(\lambda^2 - 4\lambda)$) :

$$
\lambda_1 = 4, \quad \lambda_2 = 4, \quad \lambda_3 = 0
$$

On observe ici une valeur propre répétée $\lambda = 4$ avec une multiplicité algébrique de $m_{\lambda=4} = 2$, et une valeur propre $\lambda = 0$ avec $m_{\lambda=0} = 1$.

Nous allons calculer les espaces propres pour vérifier si la multiplicité géométrique est égale à la multiplicité algébrique.

Pour la valeur propre répétée $\lambda = 4$ :
Nous cherchons l'espace propre $E_4$ en résolvant $(A - 4I)v = 0$ pour $v = (x, y, z)^\top$ :

$$
(A - 4I)v = \begin{pmatrix} 0 & 0 & 0 \\ 0 & -2 & -2 \\ 0 & -2 & -2 \end{pmatrix} \begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
$$

Ce système nous donne l'équation $-2y - 2z = 0 \implies y = -z$. La variable $x$ n'apparaît pas dans les contraintes, elle est donc libre, tout comme $z$.
Puisque nous avons deux variables libres ($x$ et $z$), nous pouvons former deux vecteurs de base linéairement indépendants :

$$
v_1 = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} \quad (\text{avec } x=1, z=0), \quad
v_2 = \begin{pmatrix} 0 \\ -1 \\ 1 \end{pmatrix} \quad (\text{avec } x=0, z=1)
$$

La dimension de l'espace propre $E_4$ est donc 2.

$$
\text{Multiplicité Géométrique } (g_{\lambda=4}) = \dim(E_4) = 2.
$$

On constate que $g_{\lambda=4} = m_{\lambda=4}$.

Pour la valeur propre simple $\lambda = 0$ :
Nous cherchons l'espace propre $E_0$ en résolvant $(A - 0I)v = 0$ :

$$
Av = \begin{pmatrix} 4 & 0 & 0 \\ 0 & 2 & -2 \\ 0 & -2 & 2 \end{pmatrix} \begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
$$

Ce système impose $4x = 0 \implies x = 0$. La deuxième ligne donne $2y - 2z = 0 \implies y = z$. En posant $z=1$, on obtient $y=1$.

$$
v_3 = \begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix}
$$

La dimension de l'espace propre $E_0$ est 1.

$$
\text {Multiplicité Géométrique}  (g_{\lambda=0}) = \dim(E_0) = 1.
$$

Ici aussi, $g_{\lambda=0} = m_{\lambda=0}$.

Pour chaque valeur propre, la multiplicité géométrique est égale à la multiplicité algébrique ($g_\lambda = m_\lambda$). Cela signifie que malgré la présence d'une racine répétée, la matrice possède un nombre suffisant de vecteurs propres linéairement indépendants ($v_1, v_2, v_3$).

```{admonition} Lemme (Rang du produit de matrices)
:class: note
Soient $A \in \mathbb{R}^{m \times n}$ et $B \in \mathbb{R}^{n \times k}$. Le rang de leur produit $AB \in \mathbb{R}^{m \times k}$ satisfait l'inégalité suivante :

$$
\operatorname{rang}(AB) \leq \min\{\operatorname{rang}(A), \operatorname{rang}(B)\}.
$$
```

```{admonition} Démonstration
Considérons le produit matriciel $AB$ :

- Chaque ligne de $AB$ est une combinaison linéaire des lignes de $B$, ce qui implique que l'espace ligne de $AB$ est contenu dans celui de $B$. Par conséquent, $\operatorname{rang}(AB) \leq \operatorname{rang}(B)$.
- De même, chaque colonne de $AB$ est une combinaison linéaire des colonnes de $A$, donc l'espace colonne de $AB$ est contenu dans celui de $A$. D'où $\operatorname{rang}(AB) \leq \operatorname{rang}(A)$.

En combinant ces observations, on conclut que $\operatorname{rang}(AB) \leq \min\{\operatorname{rang}(A), \operatorname{rang}(B)\}$.
```

```{admonition} Lemme (Rang des matrices symétriques)
:class: note
Si $A$ est une matrice réelle symétrique $n \times n$,  alors  $\operatorname{rang}(A)$ est égal au nombre total de valeurs propres non nulles de $A$. De plus, l'espace colonne $C(A)$ est le sous-espace linéaire engendré par les vecteurs propres de $A$ correspondant à ses valeurs propres non nulles.
```

```{admonition} Démonstration
Toute matrice symétrique $A$ peut être exprimée sous sa forme spectrale $A = Q \Lambda Q^\top$, où $Q$ est une matrice orthogonale et $\Lambda$ est une matrice diagonale contenant les valeurs propres de $A$. En utilisant le lemme précédent, nous procédons comme suit :

- À partir de $A = Q \Lambda Q^\top$, on a $\operatorname{rang}(A) \leq \operatorname{rang}(Q \Lambda) \leq \operatorname{rang}(\Lambda)$.
- À partir de $\Lambda = Q^\top A Q$, on a $\operatorname{rang}(\Lambda) \leq \operatorname{rang}(Q^\top A) \leq \operatorname{rang}(A)$.

Cela implique que $\operatorname{rang}(A) = \operatorname{rang}(\Lambda)$, ce qui est égal au nombre total de valeurs propres non nulles de $A$.
```

```{admonition} Exemple
:class: seealso
Soit la matrice symétrique réelle :

$$
A = \begin{bmatrix}
    8  & -2 & 2 \\
    -2 & 5  & 4 \\
    2  & 4  & 5
\end{bmatrix}
$$
```

Le polynôme caractéristique $P_A(\lambda) = -\lambda(\lambda - 9)^2$ donne les valeurs propres :

$$
\lambda_1 = 0 \quad (m_0 = 1), \quad \lambda_{2,3} = 9 \quad (m_9 = 2).
$$

Les espaces propres associés sont :

- Pour $\lambda = 0$ : $E_0 = \operatorname{vect}(v_1)$ avec $v_1 = (1, 2, -2)^\top$. ($g_0 = m_0 = 1$).
- Pour $\lambda = 9$ : $E_9 = \operatorname{vect}(v_2, v_3)$ avec $v_2 = (-2, 1, 0)^\top, v_3 = (2, 0, 1)^\top$. ($g_9 = m_9 = 2$).

La matrice est diagonalisable car les multiplicités géométriques et algébriques sont égales.

Les valeurs propres sur la diagonale de $\Lambda$ sont $\{0, 9, 9\}$. Puisqu'il y a exactement 2 valeurs non nulles, le rang est :

$$
\operatorname{rang}(A) = 2.
$$

Comme indiqué dans le lemme, l'espace colonne est engendré par les vecteurs propres correspondant aux valeurs propres non nulles.
Ici, le vecteur $v_1$ correspond à $\lambda=0$, donc il ne contribue pas à l'espace colonne.
L'espace colonne est donc engendré uniquement par $v_2$ et $v_3$ :

$$
C(A) = \operatorname{vect}(v_2, v_3).
$$

Cet exemple montre que, même en présence de valeurs propres répétées, le rang de la matrice symétrique est entièrement déterminé par les sous-espaces propres associés aux valeurs propres non nulles.

*Source : Department of Mathematics, University of Texas at Austin, "Singular Value Decomposition" — voir [Références](references.md).*

```{admonition} Remarque (Unicité de la Théorème Spectrale)
:class: hint
Il est important de noter que la Théorème Spectral d'une matrice n'est généralement pas unique. L'une des raisons est la présence de valeurs propres répétées.

Lorsque deux ou plusieurs valeurs propres $\lambda_i$ et $\lambda_j$ sont identiques, l'échange de leurs vecteurs propres correspondants dans la matrice orthogonale $Q$ résulte en une décomposition différente qui reste mathématiquement valide et équivalente.
Cependant, les *sous-espaces propres* associés à chaque valeur propre restent fixes. Cela signifie que, bien que le choix des vecteurs propres au sein de chaque sous-espace puisse varier, la décomposition en termes de sous-espaces propres est unique.

Dans l'exercice que nous venons de résoudre avec la matrice :

$$
A = \begin{bmatrix}
    8  & -2 & 2 \\
    -2 & 5  & 4 \\
    2  & 4  & 5
\end{bmatrix}
$$

nous avons observé une valeur propre répétée $\lambda = 9$ de multiplicité 2.

- **Non-unicité de la base :** Nous avons choisi une base orthonormée $\{u_2, u_3\}$ pour construire la matrice $Q$. Cependant, toute autre paire de vecteurs orthonormés appartenant à ce plan serait également un choix valide. Par exemple, nous pourrions effectuer une rotation des vecteurs $u_2$ et $u_3$ dans ce plan. Ceci produirait une matrice orthogonale différente $Q'$, mais l'égalité $A = Q'\Lambda Q'^\top$ resterait vérifiée.
- **Unicité du sous-espace :** Bien que l'on puisse modifier les vecteurs $u_2$ et $u_3$, le sous-espace propre $E_9$ (le plan lui-même) est unique et invariant.

En conclusion, pour les matrices symétriques réelles, bien que la matrice $Q$ ne soit pas unique en présence de valeurs propres multiples, la structure géométrique des sous-espaces propres est, quant à elle, parfaitement déterminée.
```
