= Bessel's Correction Term

#let avg(seq, len) = {
  $overline(attach(#seq, tl: #len))$
}
#let Chebyshev = {
  link("https://en.wikipedia.org/wiki/Chebyshev%27s_inequality#Probabilistic_statement")[Chebyshev's Inequality]
}

Let ${ X_i }$ be i.i.d. Denote expected value of random variable $X$ by $EE bracket.l X bracket.r$ and its variance by $VV bracket.l X bracket.r$. Denote generic average $#avg[$alpha$][$k$] := 1/k limits(sum)_(i=1)^k alpha_i$ and let $Y_i := #avg[$X$][$i$]$.

A naive way of estimating variance is calculating $1/k limits(sum)_(i=1)^k (x_i - #avg[x][k])^2$ and hoping it being close to $VV bracket.l X bracket.r = EE bracket.l (X - mu_X)^2 bracket.r$. Formally speaking, we're indeed saying "I want the _joint_ random variable $Z := limits(sum)_(i=1)^k (X_i - Y_k)^2$ to be s.t. its expected value $EE bracket.l Z bracket.r = EE bracket.l limits(sum)_(i=1)^k (X_i - Y_k)^2 bracket.r$ _equal_ to $k dot VV bracket.l X bracket.r$", _plus_ that "the oneshot sample of joint random variable $Z$ being close to its expected value". The latter of which is bounded by e.g. #Chebyshev, which describes the mass away from the expected value.

For the former, though, indeed we have:

$
  EE bracket.l limits(sum)_(i=1)^k (X_i - Y_k)^2 bracket.r
  &= EE bracket.l limits(sum)_(i=1)^k (1/k (limits(sum)_(j=1)^k X_i) - Y_k)^2 bracket.r \
  &= EE bracket.l limits(sum)_(i=1)^k (1/k (limits(sum)_(j=1)^k X_i) - 1/k (limits(sum)_(j=1)^(k) X_j))^2 bracket.r \
  &= 1/(k^2) EE bracket.l limits(sum)_i (limits(sum)_j (X_i - X_j))^2 bracket.r \
  &= 1/(k^2) EE bracket.l limits(sum)_i (k X_i - (limits(sum)_(j=1)^k X_j))^2 bracket.r \
  &=^"i.i.d." 1/k EE bracket.l (k X_1 - (limits(sum)_(j=1)^k X_j))^2 bracket.r \
  &= 1/k EE bracket.l (k X_1 - X_1 - X_2 - dots - X_k)^2 bracket.r \
  &= 1/k EE mat(
    align: #left,
    gap: #0.5em,
    delim: "[",
    k^2 X_1^2, -k X_1 X_1, -k X_1 X_2, -dots, - k X_1 X_k;
    -k X_1 X_1, +X_1^2, +X_1 X_2, +dots, +X_1 X_k;
    -k X_1 X_2, +X_1 X_2, +X_2 X_2, +dots, +X_2 X_k;
    dots.v;
    -k X_1 X_k, +X_1 X_k, +X_2 X_k, +dots, +X_k X_k;
  ) \
  &=^"i.i.d." 1/k ((k^2+k-2k) EE bracket.l X^2 bracket.r + (k^2 - k) EE bracket.l X X^prime bracket.r - 2 k (k-1) EE bracket.l X X^prime bracket.r)
$

The first term is collecting those of the same random variables: the main diagonoal, one on the first row $-k X_1 X_1$, and its mirror copy on the first column.

The second term is collecting all the terms in the cofactor matrix expanded at $(1, 1)$ sans its main diagonal: the cofactor matrix being $k times k$, we got $k^2 - k$ terms; they are i.i.d. so we may dub 'em simply as $EE bracket.l X X^prime bracket.r$.

Finally there are some terms on the first row and the first column: each term being $-k dot EE bracket.l X X^prime bracket.r$, and there are $(k-1)$ of them both on the first column and on the first row.

$
  & = (k^2-k)/k (EE bracket.l X^2 bracket.r - EE bracket.l X X^prime bracket.r) \
  & =^"i.i.d." (k-1) (EE bracket.l X^2 bracket.r - EE bracket.l X bracket.r EE bracket.l X^prime bracket.r) \
  & =^"i.i.d." (k-1) (EE bracket.l X^2 bracket.r - EE bracket.l X bracket.r^2) \
  & = (k-1) (VV bracket.l X bracket.r)
$

Meaning $VV bracket.l X bracket.r = 1/(k-1) EE bracket.l limits(sum)_(i=1)^k (X_i - Y_k)^2 bracket.r$. In other words, to state "huh we got samples of some i.i.d. ${ X_i }$, we should be able to estimate $VV bracket.l X bracket.r$", we're _assuming_ oneshot sample of some _joint_ (say, $Z := limits(sum)_(i=1)^k (X_i - Y_k)^2$) being close to the corresponding expected value (i.e. that mass being close to expected value, via tools like #Chebyshev), and it's this assumption and formal manipulation of expected values that brings us to Bessel's correction term.
