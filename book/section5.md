# 5. Approximation Matricielle

La Décomposition en Valeurs Singulières permet d'exprimer toute matrice $A$ comme une somme de matrices de rang 1, chacune étant pondérée par sa propre valeur singulière ($\sigma_i$). L'ordre décroissant des valeurs singulières détermine directement le niveau de dominance de ces composantes sur la matrice.

Grâce à cette structure, il est possible d'approximer efficacement la matrice $A$ en ne conservant qu'un nombre limité de termes (les $k$ premiers) de cette somme. Cette approche nous conduit naturellement au concept fondamental d'approximation matricielle de rang faible (*low-rank approximation*), qui permet de préserver la grande majorité de l'information structurelle de la matrice d'origine.

Pour une matrice $A$ de ${rang}(A) = r$, l'approximation de rang $k$ (où $k < r$) est définie par l'expression suivante :

\begin{equation}
    \widehat{A}(k) = \sum_{i=1}^{k} \sigma_i u_i v_i^{\top}
\end{equation}

```{admonition} Théorème (Eckart–Young Théorème - Norme Spectral)
:class: important
Soit $A \in \R^{m \times n}$ une matrice de rang $r$, et soit $B \in \R^{m \times n}$ une matrice de rang $k$, l'inégalité suivante est toujours satisfaite :

\begin{equation}
    \|A - B\|_2 \geq \|A - A_k\|_2 = \sigma_{k+1}
\end{equation}

Ce théorème affirme que la matrice $A_k$, obtenue par la DVS, est la meilleure approximation de rang $k$ de la matrice $A$ au sens de la norme spectrale.
```

```{admonition} Démonstration
Supposons par l'absurde qu'il existe une matrice $B$ de rang $k$ telle que :

$$
\|A - B\|_2 < \|A - A_k\|_2 = \sigma_{k+1}
$$

Pour tout vecteur non nul $w \in \R^n$, nous aurions $\|(A - B)w\|_2 < \sigma_{k+1}\|w\|_2$. En particulier, pour tout $w \in \ker(B)$ (où $Bw = 0$), cela impliquerait :

$$
\|Aw\|_2 < \sigma_{k+1}\|w\|_2
$$ (eq:1)

D'autre part, considérons le sous-espace $V_{k+1} = \text{span}\{v_1, \dots, v_{k+1}\}$. Pour tout $w \in V_{k+1}$, $w$ peut s'écrire comme $w = \sum_{i=1}^{k+1} c_i v_i$. Par l'orthogonalité des vecteurs $v_i$ et puisque $\sigma_1 \ge \dots \ge \sigma_{k+1}$, nous avons :

$$
\|Aw\|_2 = \left\| \sum_{i=1}^{k+1} c_i \sigma_i u_i \right\|_2 = \sqrt{\sum_{i=1}^{k+1} c_i^2 \sigma_i^2} \ge \sigma_{k+1} \sqrt{\sum_{i=1}^{k+1} c_i^2} = \sigma_{k+1}\|w\|_2
$$ (eq:2)

D'après le théorème du rang, $\dim(\ker(B)) = n - k$ et nous savons que $\dim(V_{k+1}) = k + 1$. La somme de leurs dimensions étant $(n-k) + (k+1) = n+1 > n$, il existe nécessairement un vecteur $w \neq 0$ tel que $w \in \ker(B) \cap V_{k+1}$.

Cependant, pour ce vecteur $w$, les inégalités {eq}`eq:1` et {eq}`eq:2` sont contradictoires. Cette contradiction prouve qu'aucune matrice $B$ de rang $k$ ne peut surpasser l'approximation $A_k$.
```

```{admonition} Théorème (Eckart–Young Théorème - Norme de Frobenius)
:class: important
Soit $A \in \mathbf{R}^{m \times n}$ une matrice de rang $r$, et soit $B \in \mathbf{R}^{m \times n}$ une matrice de rang $k$, l'inégalité suivante est toujours satisfaite :

$$
\|A - B\|_F^2 \geq \|A - A_k\|_F^2 = \sum_{i=k+1}^{r} \sigma_i^2
$$

Ce théorème affirme que la matrice $A_k$, obtenue par la DVS, est la meilleure approximation de rang $k$ de la matrice $A$ au sens de la norme de Frobenius.
```

```{admonition} Démonstration
Toute matrice $B \in \mathbf{R}^{m \times n}$ de rang $k$ peut s'écrire sous la forme d'une factorisation $k$ $B = CR$, où $C \in \mathbf{R}^{m \times k}$ et $R \in \mathbf{R}^{k \times n}$.

On peut toujours paramétrer une telle factorisation en introduisant une décomposition DVS de $B$ :
$B = U_B \Sigma_B V_B^{\top}$, et en posant
$C = U_B \Sigma_B$ et $R = V_B^{\top}$.
Dans cette écriture, on a bien $RR^{\top} = I_k$ et les colonnes de $C$ sont orthogonales, avec
$C^{\top}C = \Sigma_B^2$ diagonale.

En utilisant l'identité $\|X\|_F^2 = \text{tr}(X X^{\top})$, nous pouvons développer la fonction d'erreur $E = \|A - CR\|_F^2$ sous la forme suivante :

\begin{align}
    E & = \text{tr}\left( (CR - A)(CR - A)^{\top} \right) \nonumber                                                           \\
      & = \text{tr}(C R R^{\top} C^{\top}) - \text{tr}(C R A^{\top}) - \text{tr}(A R^{\top} C^{\top}) + \text{tr}(A A^{\top})
\end{align}

Afin d'analyser la structure des solutions de rang $k$, on considère les conditions stationnaires par rapport à $C$ et $R$.

**Remarque 1.**
*Pour le calcul des dérivées matricielles impliquant l'opérateur de trace ($\text{tr}$), nous utilisons les identités standards suivantes :*

\begin{align*}
    \frac{\partial}{\partial X} \text{tr}(XA)         & = A^{\top}       \\
    \frac{\partial}{\partial X} \text{tr}(XAX^{\top}) & = XA^{\top} + XA \\
    \frac{\partial}{\partial X} \text{tr}(AX^{\top})  & = A
\end{align*}

En dérivant par rapport à $C$, on obtient :

$$
\frac{\partial E}{\partial C} = 2C R R^{\top} - 2A R^{\top}
$$

De même, la dérivation par rapport à $R$ donne :

$$
\frac{\partial E}{\partial R} = 2C^{\top}C R - 2C^{\top}A
$$

Aux points stationnaires, on a donc :

$$
C = A R^{\top}, \qquad C^{\top}A = (C^{\top}C)R
$$

En combinant ces deux relations, on obtient :

$$
(A^{\top}A)R^{\top} = R^{\top}(C^{\top}C)
$$

Ce qui montre que les colonnes de $R^{\top}$ engendrent un sous-espace invariant de $A^{\top}A$, donc elles coïncident avec les vecteurs singuliers à droite associés aux $k$ plus grandes valeurs singulières, i.e. $R^{\top} = V_k^{\top}$.

Réciproquement, on obtient :

$$
(AA^{\top})C = C(C^{\top}C)
$$

Ce qui implique que les colonnes de $C$ appartiennent au sous-espace engendré par les vecteurs singuliers à gauche $U_k$ de $A$.

Ainsi, toute solution optimale de rang $k$ est de la forme :

\begin{equation}
    B = CR = U_k \Sigma_k V_k^{\top} = \sum_{i=1}^{k} \sigma_i u_i v_i^{\top}.
\end{equation}

L'erreur associée s'écrit alors :

\begin{equation}
    \|A - B\|_F^2 = \sum_{i \notin S} \sigma_i^2.
\end{equation}

Pour minimiser cette quantité, il est nécessaire de choisir les indices correspondant aux plus grandes valeurs singulières de $A$. Ainsi, on prend

\begin{equation}
    S = \{1, 2, \dots, k\},
\end{equation}

ce qui donne la meilleure approximation.

On obtient alors :

\begin{equation}
    \|A - B\|_F^2 = \sum_{i=k+1}^{r} \sigma_i^2, \quad \text{et} \quad B = A_k = U_k \Sigma_k V_k^{\top}.
\end{equation}

Cela conclut la démonstration du théorème d'Eckart–Young. $\square$
```
