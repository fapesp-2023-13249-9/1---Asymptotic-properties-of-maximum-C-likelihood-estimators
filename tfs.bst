\documentclass[]{interact}

\usepackage{epstopdf}% To incorporate .eps illustrations using PDFLaTeX, etc.
\usepackage{subfigure}% Support for small, `sub' figures and tables
\usepackage{multirow}
\usepackage{float}
\usepackage[numbers,sort&compress]{natbib}% Citation support using natbib.sty
\bibpunct[, ]{[}{]}{,}{n}{,}{,}% Citation support using natbib.sty
\renewcommand\bibfont{\fontsize{10}{12}\selectfont}% Bibliography support using natbib.sty
\newcommand{\f}{\operatorname}


\theoremstyle{plain}% Theorem-like structures
\newtheorem{theorem}{Theorem}[section]

\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{corollary}[theorem]{Corollary}
\newtheorem{proposition}[theorem]{Proposition}

\theoremstyle{definition}
\newtheorem{definition}[theorem]{Definition}
\newtheorem{example}[theorem]{Example}

\theoremstyle{remark}
\newtheorem{remark}{Remark}
\newtheorem{notation}{Notation}
\newcommand{\R}{\mathbb{R}}
\begin{document}

\articletype{ARTICLE TEMPLATE}

\title{A modified closed-form maximum likelihood estimator}

\author{ Pedro Luiz Ramos$^{\rm a}$$^{\ast}$\thanks{$^\ast$Corresponding author. Email: pedrolramos@usp.br
\vspace{6pt}}, Eduardo Ramos, Francisco Louzada and Francisco Rodrigues$^{\rm a}$ \\\vspace{6pt}  $^{a}${Institute of Mathematical Science and Computing, University of S\~ao Paulo, S\~ao Carlos, Brazil} }


\maketitle

\begin{abstract}
In progress\end{abstract}

\begin{keywords}
Keywords 1; Keywords 2; Keywords 3; Keywords 4.
\end{keywords}

\section{Introduction}


Statistics is widely concerned in proposing and improving inferential methods. Introduced by Ronald Fisher, the maximum likelihood method is one of the most well-know and used inferential procedure to estimate the unknown parameters of given model. 
Alternative methods to the MLE have been considered or proposed in the literature, some are based on moments \cite{hosking1990moments}, percentile \cite{kao1,kao2}, product of spacings \cite{Cheng2}, or goodness of fit measures, to list a few. Although, alternative frequentist methods are popular nowadays, the MLEs is the most widelly used due to its flexibility in include additional complexity (such as random effects, covariates, censoring) as well as its good properties such as asymptotically efficiency, consistent and invariance under one-to-one transformation. These properties are achieve when the MLEs satisfies some regularity conditions (include reference).


It is well known that the maximum likelihood (ML) estimators does not always return closed-form expressions for most of the common models. In these cases, numerical methods such as Newton-Rapshon or its variants are usually considered to find the values that maximizes the likelihood function. Important variants of the maximum likelihood estimator such as profile \cite{murphy2000profile}, conditional \cite{andersen1970asymptotic}, and marginal likelihoods \cite{cox1975partial} have been presented in order to  to eliminate nuisance parameters and decrease the computational cost. Another important procedure to achieve the MLEs is the expectation–maximization (EM) algorithm which involves the addition of unobserved latent variables jointly with unknown parameters. The expectation and maximization steps also involves, in most cases, the use of numerical methods that may have an significant computational time. Although, in the recent decade we have overcome computational limitation in terms of software and hardware. On the other hand, in many situations there is an need of using closed-form estimators to estimate the unknown parameters, for instance, dealing with embed technology when small components need to compute the estimates without the use of implemented maximization procedures, or cases when dealing with large amount of data that are used in real time decision and need to be computed quickly.

In this paper, we discuss a modification in the maximum likelihood method that, under certain conditions, allows us to obtain closed-form expressions for the parameters of the model. The proposed estimator has important properties such as consistency and asymptotic normality. A formal proof is provided where the conditions do not need to be verified if the distribution belongs to the exponential family. Overall, to find the closed-form estimator, the method depends on a generalized version of the standard distribution, and based on its likelihood equations we obtain the estimators. The new method is illustrated in the Gamma, Beta, and Nakagami distribution, where the standard ML does not have closed-form expressions, and numerical methods or approximations are necessary to find the solutions for the cited models. Hence our approach does not require iterative numerical methods and computationally the work required with our estimators is less arduous than that required in ML estimators. The remainder of this paper is organized as follows. Section 2 presents the new modified maximum likelihood estimator and its properties. Section 3 considers the application in the gamma, Nakagami, and Beta distributions. Finally, Section 6 summarizes the study.




\section{Modified Maximum Likelihood Estimator}

The new suggested method can be applied to obtain closed-form expressions if some conditions are fulfilled. It is describe to a distribution $f(t\ ;\ \theta)$  and cumulative density function (c.d.f.) $F(t\ ;\ \theta)$.

Let $\mathcal{A}\subset \R^{s}$ and $\mathcal{B}\subset \R^{r}$ be open sets with $s>0$ and $0\leq r\leq s$, suppose that $f(t\ ;\ \theta,\alpha): E\to \R$ is a distribution in $t$ for each fixed $(\theta,\alpha)\in \mathcal{A}\times \mathcal{B}$, where $E\subset \R$ is an open set.

Then, letting $l_n(\theta,\alpha\ ;\ t)=\sum_{i=1}^n \log\left(f(t_i\ ;\ \theta,\alpha)\right)$ for all $n>0$ and $t\in E^n$ and given a fixed $\alpha_0\in \mathcal{B}$, we consider the hybrid likelihood equations over $\theta$ given by
\begin{equation}\label{hybrid}
\begin{aligned}
\frac{\partial}{\partial \theta_i} l_n\left(\theta,\alpha_0\ ;\ X\right) = 0,\ 1\leq j\leq s-r\\
\frac{\partial}{\partial \alpha_j} l_n\left(\theta,\alpha_0\ ;\ X\right)= 0,\ 1\leq j\leq r.
\end{aligned}
\end{equation}

Given $X_1,X_2,\cdots$ iid, with a  distribution $f(t\ ;\ \theta_0,\alpha_0)$, our objective shall be that of giving conditions such that there exists $\hat{\theta}_n(X)$ satisfying the hybrid likelihood equations above such that $\hat{\theta}_n(X)$ is a consistent estimator for $\theta_0$.

The likelihood equations are clearly an special case of the above equations when $r=0$.

Now, in order to formulate our theorem, given a fixed $\alpha_0\in \mathcal{A}$, define $I(\theta)=\left(I_{j,k}(\theta)\right)\in M_s(\R)$ as usual by
 \begin{equation*}
 \begin{aligned}
 I_{i,j}(\theta)=cov_{\theta_0} \left[\frac{\partial\,  \log\, f(X\; ;\; \theta,\alpha_0)}{\partial \theta_i},\frac{\partial\,  \log\, f(X\; ;\; \theta,\alpha_0)}{ \partial \theta_j}\right]
 \end{aligned}
 \end{equation*}
for all $1\leq i\leq s$, $1\leq j\leq s$ and define $J(\theta)=\left(J_{i,j}(\theta)\right)\in M_{s}(\R)$, by
 \begin{equation}\label{eqj}
 \begin{aligned}J_{i,j}(\theta)=
 cov_{\theta_0} \left[\frac{\partial\,  \log\, f(X\; ;\; \theta,\alpha_0)}{\partial \theta_i},\frac{\partial\,  \log\, f(X\; ;\; \theta,\alpha_0)}{ \partial \beta_j}\right]
 \end{aligned}
 \end{equation}
 for all $1\leq i\leq s$, $1\leq j\leq s$, where $(\beta_1,\cdots,\beta_s)=(\theta_1,\cdots,\theta_{s-r},\alpha_1,\cdots,\alpha_r)$.
\begin{theorem}\label{princ} Let $X=\left(X_1, X_2, \cdots, X_n\right)$ be iid with density $f(t\ ;\ \theta_0,\alpha_0): E\to \R$ and suppose $(\hat{\theta}_{1n}(X),\cdots,\hat{\theta}_{sn}(X))$ satisfy the following properties for all $n\geq 2$.
\begin{itemize}
\item[(A)] With probability $1$ over $w$, the hybrid-likelihood equations has  $\hat{\theta}(X)$ as its unique solution.
\item[(B)]  $J(\theta_0)$, as defined in (\ref{eqj}), is invertible.
\item[(C)] $f(t\ ;\ \phi)$ is a $C^2$ function in $\phi=(\phi_1,\cdots,\phi_{s+r})=(\theta_1,\cdots,\theta_s,\alpha_1,\cdots,\alpha_r)$ for all $t\in E$ and moreover:
\begin{equation*}\int_{E}\frac{\partial
 f(t\ ;\ \phi)}{\partial \phi_i}dt =0\mbox{ and }\int_{E}\frac{\partial^2
 f(t\ ;\ \phi)}{\partial \phi_i\partial \phi_j}dt =0
\end{equation*}

%\begin{equation*}
% \begin{aligned}
% & E_{\phi}\left[\frac{\partial\,  \log\, f(X\; ;\; \phi)}{\partial \phi_i}\right]=0\mbox{ and }\\
% & E_{\phi}\left[\frac{\partial\,  \log\, f(X\; ;\; \phi)}{\partial \phi_i}\cdot \frac{\partial\,  \log\, f(X\; ;\; \phi)}{\partial \phi_j}\right]=E_{\phi}\left[-\frac{\partial^2\,  \log\, f(X\; ;\; \phi)}{\partial \phi_i\partial \phi_j}\right]
% \end{aligned}
% \end{equation*}
for $1\leq i\leq s+r$ and $1\leq j\leq s+r$.
\item[(D)] There exist measurable functions $M_{ij}$ and an open set $\Theta_0\subset \Theta$ containing  $\phi_0=(\theta_0,\alpha_0)$ such that
\begin{equation*}
 \begin{aligned}
 \left|\frac{\partial^2\,  \log\, f(X\; ;\; \phi)}{\partial \phi_i \partial \phi_j}\right|\leq M_{ij}(X)\mbox{ and }E_{\phi_0}\left[M_{ij}(X)\right]<\infty
 \end{aligned}
 \end{equation*}
for all $\phi\in \Theta_0$, $1\leq i\leq s+r$ and $1\leq j\leq s+r$.
\end{itemize}
Then:
\begin{itemize}
    \item[I)] $\hat{\theta}_n(X)$ is a strongly consistent estimator for $\theta_0$.
    \item[II)]
$\sqrt{n}(\hat{\theta}_n-\theta_0)\overset{D}{\to} N_s\left(0,(J(\theta_0))^{-1}\right)$.
\end{itemize}
\end{theorem}

Usually, for $n\geq 2$, the hybrid likelihood equations have a unique solution unless $X_1=X_2=\cdots=X_n$. Since in this paper the parameter space is continuous it easy to see that such equality has probability zero of occurring, and thus condition $(A)$ is satisfied under such conditions. Moreover, we notice that condition $(B)$ is a simple generalization of a usual condition on the literature which asks that the Fischer information matrix $I(\theta_0)$ be invertible since for $r=0$ it follows that $J(\theta_0)=I(\theta_0)$.

Condition (C) is equivalent to asking that the probability density equation
\begin{equation}\label{eqdens}\int_{E}
 f(t\ ;\ \theta,\alpha)dt =1\mbox{ for all }(\theta,\alpha)\in \Theta\times \mathcal{A}
 \end{equation}
 can be differentiated twice under the integral sign over the parameters $\theta$ and $\alpha$. We notice moreover that this condition is a standard on the literature, and is satisfied for a large set of distributions. Finally, we notice that condition $(D)$ is also standard in the literature and is satisfied for a large set of distributions (see Lehman).
 
\begin{corollary}\label{princ} Let $X=\left(X_1, X_2, \cdots, X_n\right)$ be iid with density $f(t\ ;\ \theta_0,\alpha_0): E\to \R$ where
\begin{equation*}f(t ; \theta,\alpha)=h(\theta,\alpha)\exp(-\phi(\theta,\alpha)T(X))
\end{equation*}
and suppose $\hat{\theta}(X)=(\hat{\theta}_{1n}(X),\cdots,\hat{\theta}_{sn}(X))$ satisfy the following properties for all $n\geq 2$.
\begin{itemize}
\item[(A)] With probability $1$ over $w$, the hybrid-likelihood equations has  $\hat{\theta}(X)$ as a solution.
\item[(B)]  $J(\theta_0)$, as defined in (\ref{eqj}), is invertible.
\end{itemize}

Then the following hold:
\begin{itemize}
    \item[I)] $\hat{\theta}_n(X)$ is a strongly consistent estimator for $\theta_0$.
    \item[II)]
$\sqrt{n}(\hat{\theta}_n-\theta_0)\overset{D}{\to} N_s\left(0,(J(\theta_0))^{-1}\right)$.
\end{itemize}
\end{corollary}

\section{Examples}

\vspace{0.3cm}

We illustrate the discussion of the previous section by applying the proposed approach in Gamma, Nakagami-m and Beta distribution. Since these are all part of the exponential family, following Lehman, it follows that conditions $(C)$ and $(D)$ are automatically satisfied for these models. The examples are for illustration, so we shall not present their backgrounds. The standard ML for the cited distributions are widely discussed in statistical distributions which shows that no closed-form expression can be achieved using the ML method.

\noindent\textbf{Example 1:} Let us consider that $X$  is a random variable (RV) with gamma distributions and probability density function (PDF) given by
\begin{equation}\label{fdpgamma}
f(t|\phi,\lambda)=\frac{1}{\Gamma(\phi)}\left(\frac{\phi}{\lambda} \right)^\phi t^{\phi-1}\exp\left(\frac{\phi}{\lambda} t \right),
\end{equation}
where $\phi>0$ is the shape parameter,  $\lambda >0$ is the scale parameter and $\Gamma(\alpha)=\int_{0}^{\infty}{e^{-x}x^{\alpha-1}dx}$ is the  gamma function.
The modified maximum likelihood estimators are given by
\begin{equation}\label{closmdfy}  
\hat\phi=\cfrac{n\sum_{i=1}^n X_i}{\left(n\sum_{i=1}^{n}X_i \log\left(X_i\right) - \sum_{i=1}^n X_i \sum_{i=1}^n \log\left(X_i\right) \right)} 
\end{equation}
and
\begin{equation}\label{closmdfy2}
\hat\lambda=\frac{1}{n}\sum_{i=1}^n{X_i}.
\end{equation}


This estimator can be obtained by considering that the generalized model $\mathcal{M}_\alpha$ follows a generalized gamma distribution with likelihood function given by
\begin{equation}\label{verognk1}
L(\phi,\lambda,\alpha;\boldsymbol{t})=\frac{\alpha^n}{\Gamma(\phi)^n}\left(\frac{\phi}{\lambda} \right)^{n\phi}\left\{\prod_{i=1}^n{t_i^{\alpha\phi-1}}\right\}\exp\left(-\frac{\phi}{\lambda} \sum_{i=1}^n t_i^\alpha \right). \end{equation}

The maximum likelihood estimates of the parameters are obtained by solving the following likelihood equations 
\begin{equation}\label{verogg21}  
\hat\phi=\frac{n}{\left(\frac{1}{\hat\lambda}\sum_{i=1}^n t_i^\alpha \log\left(t_i^\alpha\right) -\sum_{i=1}^n \log\left(t_i^\alpha\right) \right)} \ , 
\end{equation}
\begin{equation}\label{verogg22}  
\hat\lambda=\frac{1}{n}\sum_{i=1}^n{t_i^{\hat\alpha}}
\end{equation}
and the MLE for $\alpha$ will be obtained solving the non-linear equation
\begin{equation}\label{verogg23} 
\log(\hat\phi)-\psi(\hat\phi)=\log(\hat\lambda)-\frac{1}{n}\sum_{i=1}^n{\log\left(t_i^{\hat\alpha}\right)}. 
 \end{equation}

Note that the generalized model $\mathcal{M}_\alpha$ satisfies Theorem \ref{princ}, i.e., $\lambda=f_1(\boldsymbol{t},\alpha)$; $\phi=f_2(\boldsymbol{t};\lambda,\alpha)$, where $f_1(\cdot)$ and $f_2(\cdot)$ have closed-form expressions. For $\alpha=1$, we have that the Gamma distributions has closed-form estimators (see Louzada et al. \cite{louzada2019note}) which is given by
\begin{equation*}
\begin{aligned}
&\hat\lambda=\frac{1}{n}\sum_{i=1}^n{t_i} \ \ \mbox{ and }
\hat\phi=\cfrac{\sum_{i=1}^n t_i}{\left(\sum_{i=1}^{n}t_i \log\left(t_i\right) - \frac{1}{n}\sum_{i=1}^n t_i \sum_{i=1}^n \log\left(t_i\right) \right)} \, \cdot
\end{aligned}
\end{equation*}


\vspace{0.3cm}

\noindent\textbf{Example 2:} Let $T$ be a random variable that follows a Nakagami-m distributions with PDF given by
\begin{equation*}\label{fdpnk}
f(t|\mu,\mathcal{A})=\frac{2}{\Gamma(\mu)}\left(\frac{\mu}{\mathcal{A}} \right)^\mu t^{2\mu-1}\exp\left(\frac{\mu}{\mathcal{A}} t^2 \right) 
\end{equation*}
for all $t>0$, where $\mu\geq 0.5$ and $\mathcal{A}>0$.


From Example 1, the generalized model $\mathcal{M}_\alpha$ also satisfies the Proposition \ref{proposicao1}, i.e., $\lambda=f_1(\boldsymbol{t},\alpha)$; $\phi=f_2(\boldsymbol{t};\lambda,\alpha)$, where $f_1(\cdot)$ and $f_2(\cdot)$ have closed-form expressions. For $\alpha=2$ (see Ramos et al. \cite{ramos2016efficient}) the modified maximum likelihood estimator is given by
\begin{equation*}\label{vinkom}
\begin{aligned}
&\hat\lambda=\frac{1}{n}\sum_{i=1}^n{t_i^2} \ \ \mbox{ and }
\hat\phi=\cfrac{\sum_{i=1}^n t_i^{2}}{\left(\sum_{i=1}^{n}t_i^{2} \log\left(t_i^2\right) - \frac{1}{n}\sum_{i=1}^n t_i^{2} \sum_{i=1}^n \log\left(t_i^2\right) \right)} \, \cdot
\end{aligned}
\end{equation*}

Note that the approach given above can be considered for other particular cases, for instance, the Wilson-Hilferty distributions is obtained when $\alpha=3$. Hence, we can obtain closed-form estimators for cited distribution as well.

\vspace{0.3cm}

\noindent\textbf{Example 3:} Let us consider the beta distribution where the PDF is given by
\begin{equation}\label{fdpbeta}
f(x)=\frac{(\phi x)^{a-1}(1-x)^{b-1}}{\f{B}(a,b)}, \quad  0 <x<1
\end{equation}
where $\f{B}(\alpha,\phi)=\frac{\Gamma(\alpha)\Gamma(\phi)}{\Gamma(\alpha+\phi)}$ $\Gamma(y)=\int_{0}^\infty x^{y-1}e^{-y x}dx$, $\alpha>0$, and $\phi>0$.


It is well-known that both ML estimators does not have closed-form expressions and depend on transcendental functions. To obtain the modified MLE let us consider the McDonald and Richards \cite{mcdonald1987some} beta power distribution with four parameters and PDF given by
\begin{equation}\label{fdpbeta}
f(x)=\frac{\alpha\phi(\phi x)^{a\alpha-1}(1-(\phi x)^\alpha)^{b-1}}{\f{B}(a,b)}, \quad  0 <x<\frac{1}{\phi}
\end{equation}
where $a,b,\alpha,\phi>0$. The generalized form reduces to standard Beta distribution when $\alpha=\phi=1$.  The likelihood function from (\ref{fdpbeta}) is given by
\begin{equation*}\label{verognk1}
L(\boldsymbol\phi;\boldsymbol{t})=\frac{\alpha^n\phi^n}{\f{B}(\alpha,\phi)^n}\prod_{i=1}^n (\phi x_i)^{a\alpha-1}\prod_{i=1}^n (1-(\phi x_i)^\alpha)^{b-1}.\end{equation*} 

The maximum likelihood estimates of the parameters are obtained by solving: 
\begin{equation}\label{eqb1}
-n\psi(a)+n\psi(a+b)+\alpha\sum_{i=1}^{n}\log(\phi x_i)=0
\end{equation}
\vspace{-0.7cm}
\begin{equation}\label{eqb2}
-n\psi(b)+n\psi(a+b)+\sum_{i=1}^{n}\log\left(1-(\phi x_i)^{\alpha}\right)=0
\end{equation}
\vspace{-0.7cm}
\begin{equation}\label{eqb3}
\frac{n}{\alpha}+a\sum_{i=1}^{n}\log(\phi x_i)-(b-1)\sum_{i=1}^{n}\frac{(\phi x_i)^{\alpha}\log(\phi x_i)}{1-(\phi x_i)^{\alpha}}=0
\end{equation}
\vspace{-0.7cm}
\begin{equation}\label{eqb4}
\frac{\alpha}{\phi}\left(na -(b-1)\sum_{i=1}^{n}\frac{(\phi x_i)^{\alpha}}{1-(\phi x_i)^{\alpha}} \right)=0
\end{equation}

Note that considering $\alpha=\phi=1$ in (\ref{eqb1}) and (\ref{eqb2}) the maximum likelihood equations reduce to the standard MLE for the Beta distribution, which does not have closed-form expression. On the other hand, from (\ref{eqb3}) and (\ref{eqb4}) we have that
\begin{equation*}
n +\sum_{i=1}^{n}\log(x_i)-(b-1)\sum_{i=1}^{n}\frac{x_i\log(x_i)}{1-x_i}=0
\end{equation*}
\vspace{-0.7cm}
\begin{equation*}
na -(b-1)\sum_{i=1}^{n}\frac{x_i}{1-x_i}=0
\end{equation*}
which after some algebraic manipulations return the closed-form estimators 
\begin{equation}
b=n\alpha\left(\sum_{i=1}^{n}\frac{x_i}{1-x_i} \right)^{-1}+1
\end{equation}
\vspace{-0.7cm}
\begin{equation}
a=\frac{n\sum_{i=1}^{n}\frac{x_i}{1-x_i} }{n\sum_{i=1}^{n}\frac{x_i\log(x_i)}{1-x_i}-\sum_{i=1}^{n}\log(x_i)\sum_{i=1}^{n}\frac{x_i}{1-x_i}  }
\end{equation}

It is important to point out that the modified maximum likelihood estimators may have different forms depending on the generalized model. For instance, different generalized beta distributions can led to different estimators which will be showed in the next example.

\vspace{0.3cm}

\noindent\textbf{Example 4:} Let us assume that the chosen generalized beta distribution has the PDF given by
\begin{equation}\label{fdpbeta}
f(x)=\frac{(y-a)^{\alpha-1}(c-x)^{\phi-1}}{(c-a)^{\alpha+\phi-1}\f{B}(\alpha,\phi)} \quad  a <x<c
\end{equation}
From the equation above we can also obtain the MML estimators. To achieved the estimators consider the log-likelihood function from (\ref{fdpbeta}) given by
\begin{equation*}\label{verognk1}
\begin{aligned}
l(\boldsymbol\phi;\boldsymbol{x})&=(\alpha-1)\sum_{i=1}^{n}\log(x_i-\alpha)+(\phi-1)\sum_{i=1}^{n}\log(c-x_i)-n\log(\f{B}(\alpha,\phi))\\&-n\left(\alpha+\phi -1\right)\log(c-a)
\end{aligned}
\end{equation*}

The MLE are obtained from equating the score functions equal zero. Although we have to compute four likelihood equations, we will present only two that will be used to obtain the modified MLE, the equations are
\begin{equation*}
{\frac {\partial l(\boldsymbol\phi;\boldsymbol{x})}{\partial a}}=-(\alpha -1)\sum _{i=1}^{b}{\frac {1}{x_{i}-a}}\,+n(\alpha +\phi -1){\frac {1}{c-a}}=0
\end{equation*}
\vspace{-0.7cm}
\begin{equation*}
{\frac {\partial l(\boldsymbol\phi;\boldsymbol{x})}{\partial c}}=(\phi -1)\sum _{i=1}^{n}{\frac {1}{c-x_{i}}}\,-n(\alpha +\phi -1){\frac {1}{c-a}}=0
\end{equation*}
again considering $a=c=1$ and after some algebraic manipulations the closed-form estimators are
\begin{equation}
\phi=\left(\sum_{i=1}^{n}\frac{x_i}{1-x_i}\right)\left(\sum_{i=1}^{n}\frac{x_i}{1-x_i}-n-\frac{n^2}{\sum_{i=1}^{n}\frac{1}{x_i}-n}\right)^{-1}
\end{equation}
\vspace{-0.7cm}
\begin{equation}
\alpha=n\phi\left(\sum_{i=1}^{n}\frac{1}{x_i}-n\right)+1
\end{equation}


\section{Conclusions}
In Progress

\section*{Disclosure statement}

No potential conflict of interest was reported by the author(s)


\section*{Acknowledgements}

%The authors are thankful to the Editorial Board and to the reviewers for their valuable comments and suggestions which led to this improved version. 
The research was partially supported by CNPq, FAPESP and CAPES of Brazil.

\section*{Appendix}

In order to prove Theorem \ref{princ} we shall need the following lemma, which answers us the following question: If a sequence of functions $F_n$ satisfy $\lim_{n\to \infty} F_n(x_0)=0$, under which conditions can we guarantee the existence of zeros $x_n$ of $F_n$ such that $\lim_{
n\to \infty} x_n = x_0$? Since this question is central when trying proving strong consistency of estimators, we believe such lemma can be used to prove such property for a large set of estimators in the literature.

\begin{lemma}\label{clemma} Let $U\subset \R^m$ be open, let $F_i:U\to \mathbb{R}^m$ be differentiable for all $i\in \mathbb{N}$, let $x_0\in U$ and suppose that
\begin{equation*}\lim_{n\to \infty} F_n(x_0)=0\mbox{ and }\lim_{n\to \infty}F_n'(x)=G(x)\mbox{ uniformly on }U,
\end{equation*}
 where $G:U\to M_n(\mathbb{R})$ is continuous at $x_0$ and $G(x_0)$ is invertible. Then there exists $N>0$ and a sequence $\{x_n\}_{n\geq N}$ in $U$ such that
 \begin{equation}\label{lcond}\lim_{n\to \infty} x_n = x_0 \mbox{ and }F_n(x_n)=0
\mbox{ for all } n\geq N,
 \end{equation}
\end{lemma}

\begin{appendix} \begin{proof} Denote $\lambda=\left\|G(x_0)^{-1}\right\|>0$. Since $U$ is open and $G$ is continuous at $x_0$ there exists $\delta>0$ such that \begin{equation}B(x_0,\delta)\subset U\mbox{ and }\label{leq1}\left\|G(x)-G(x_0)\right\|< \frac{1}{4\lambda}\mbox{ for all }x\in B(x_0,\delta).
\end{equation}
Moreover, due to the hypothesis it follows that there exists $N>0$ such that $n\geq N$ implies that
\begin{equation}\label{leq2}
\left\|F_n(x_0)\right\|< \frac{\delta}{2\lambda}\mbox{ and }\left\| F_n'(x)-G(x)\right\|< \frac{1}{4\lambda}\mbox{ for all }x\in U.
\end{equation}
Notice the inequality in (\ref{leq1}) combined with the second inequality in (\ref{leq2}) and the triangle inequality implies that
\begin{equation}\label{leq3} \left\|F_n'(x)-G(x_0)\right\|< \frac{1}{2\lambda}\mbox{ for all }x\in B(x_0,\delta)\mbox{ and }n\geq N
\end{equation}
Moreover, if we denote $\epsilon_n=2\lambda \left\|F_n(x_0)\right\|$ for all $n\geq N$, from the first inequality in (\ref{leq2}) it follows that $\epsilon_n< \delta$ and from $B(x_0,\delta)\subset U$ we conclude that
\begin{equation*}
\bar{B}(x_0,\epsilon_n)\subset U\mbox{ for all }n\geq N.
\end{equation*}

Now, given $n\geq N$ and letting $L_n:U\to \R^m$ be defined by
$L_n(x)=x-G(x_0)^{-1} F_n(x)$ for all $x\in U$
we shall prove that $L_n(\bar{B}(x_0,\epsilon_n))\subset \bar{B}(x_0,\epsilon_n)$. Indeed, notice that from the chain rule $L_n$ is differentiable in $U$ with
\begin{equation*}
L_n'(x)=I-G(x_0)^{-1} F_n'(x)=G(x_0)^{-1}\left(G(x_0)-F_n'(x)\right)\mbox{ for all }x\in U.
\end{equation*}
Thus for all $x\in B(x_0,\delta)$ we have from (\ref{leq3}) that
\begin{equation*}
\begin{aligned}
\label{peq1}
\left\|L_n'(x)\right\|=\left\|G(x_0)^{-1}\left(G(x_0)-F_n'(x)\right)\right\|\leq \left\|G(x_0)^{-1}\right\|\left\|G(x_0)-F_n'(x)\right\|\leq\frac{1}{2}.
\end{aligned}
\end{equation*}
and due to the mean value inequality it follows that
\begin{equation}\label{ineq3}
\left\|L_n(x)-L_n(x_0)\right\|\leq \frac{1}{2}\left\|x-x_0\right\|\mbox{ for all }x\in \bar{B}(x_0,\epsilon_n)
\end{equation}
Moreover, notice that
\begin{equation}\label{peq2}\left\|L_n(x_0)-x_0\right\|=\left\|G(x_0)^{-1}F_n(x_0)\right\|\leq \left\|G(x_0)^{-1}\right\|\left\|F_n(x_0)\right\|= \frac{\epsilon_n}{2}.
\end{equation}
Thus, given $x\in \bar{B}(x_0,\epsilon_n)$ from the triangle inequality and the inequality (\ref{peq2}) we have
\begin{equation*}
\begin{aligned}
\left\|L_n(x)-x_0\right\|\leq \left\|L_n(x)-L_n(x_0)\right\|+\left\|L_n(x_0)-x_0\right\|\leq \frac{\epsilon_n}{2}+\frac{\epsilon_n}{2}=\epsilon_n,
\end{aligned}
\end{equation*}
that is, $L_n(x)\in \bar{B}(x_0,\epsilon_n)$, which proves that $L_n(\bar{B}(x_0,\epsilon_n)\subset \bar{B}(x_0,\epsilon_n)$.

Now, since $L_n$ is continuous and  $L_n(\bar{B}(x_0,\epsilon_n)\subset \bar{B}(x_0,\epsilon_n)$, from the Brouwer fixed point theorem we conclude $L_n$ has at least one fixed point $x_n$ in $\bar{B}(x_0,\epsilon_n)$ for all $n\geq N$, that is:
\begin{equation*}L_n(x_n)=x_n\Rightarrow G(x_0)^{-1}F_n(x_n)=0\Rightarrow F_n(x_n)=0\mbox{ for all }n\geq N
\end{equation*}
Moreover, notice that $\lim_{n\to \infty} \epsilon_n = \lim_{n\to \infty} 2\lambda \left\|F_n(x_0)\right\|=0$, and since $x_n\in \bar{B}(x_0,\epsilon_n)$ for $n\geq N$, it follows that $\left\|x_n-x_0\right\|\leq \epsilon_n$ for all $n\geq N$ and thus
\begin{equation*}\lim_{n\to \infty}\left\|x_n-x_0\right\|=0 \Rightarrow \lim_{n\to \infty}x_n = x_0
\end{equation*}
which concludes the proof.
%Now, to conclude the proof, we shall construct the sequence $\{A_n\}_{n\geq N}$ in $M_m(\R)$ satisfying (\ref{lcond}). Letting $F_{in}:U\to \R$ denote the $i$-th coordinate of $F_n$ for each $1\leq i\leq m$, that is, $F_n = \{F_{1n},\cdots,F_{mn}\}$, and letting $x_0=(x_{01},\cdots,x_{0m})$ and $x_n=\{x_n\}$ due to the mean value theorem, for every $i$ with $1\leq i\leq m$ there exists a point $y_{i,n}\in U$ in the line connecting $x_0$ and $x_n$ such that
%\begin{equation*} F_{in}(x_n) = F_{in}(x_0) + \left( F_{ni}'(y_{n,i})\right) (x_{in}-x_{i0}).
%\end{equation*}
%Thus $b_{ij}^n=\frac{\partial F_{in}(y_{i,n})}{\partial x_j}$ for all $i$ and $j$ it follows that
%\begin{equation*}F_n(x_n) = F_n(x_0) + B_n (x_n-x_0)
%\end{equation*}
%and from $F_n(x_n)=0$ for all $n\geq N$ we have
%\begin{equation*} F_n(x_0) = -B_n(x_n-x_0)
%\end{equation*}

\end{proof}


Using the above lemma we are ready to prove theorem \ref{princ}
\begin{proof}
$\hat{\theta}_n$ is a strongly consistent estimator for $\theta_0$:\vspace{0.3cm}

 By hypothesis $(A)$ there exists a set $\mathcal{A}_0$ of probability $1$ such that $\hat{\theta}(X(w))$ is the unique solution of the hybrid likelihood equations for all $w\in \mathcal{A}_0$.

Denoting $(\beta_1,\cdots,\beta_s)=(\theta_1,\cdots,\theta_{s-r},\alpha_1,\cdots,\alpha_r)$, let $F_n:\theta\times \mathcal{A}\subset \R^s\times \mathcal{A}\to \R^s$ be defined by $F_n = \left(F_{n1},\cdots,F_{ns}\right)$ where
\begin{equation*}
\begin{aligned}F_{nj}(\theta,w)=\frac{\partial}{\partial \beta_j} l\left(\theta,\alpha_0\ ;\ X(w)\right)
\end{aligned}
\end{equation*}
for all $\theta\in \Theta$, $1\leq j\leq s$. Notice due to the strong law of the large numbers and hypothesis $(C)$ that,  we have
\begin{equation*} \lim_{n\to \infty} F_{nj}(\theta_0,w) = E\left(\frac{\partial}{\partial \beta_j} l\left(\theta_0,\alpha_0\ ;\ X\right)\right)=\int_E \frac{\partial f(t\ ;\ \theta_0,\alpha_0)}{\partial \beta_j} dt = 0\mbox{ a.s. in }\mathcal{A}
\end{equation*}
for $1\leq j\leq s$,
or, in other words, there exists a set $\mathcal{A}_1$ of probability $1$ such that
\begin{equation}\lim_{n\to \infty} F_n(\theta_0,w) = 0\mbox{ for all }w\in \mathcal{A}_1
\end{equation}
On the other hand, due to condition $(D)$ it follows that $J(\theta)$ is continuous at $\theta_0$. Moreover, from the uniform strong law of the large numbers it follows that
\begin{equation*}\lim_{n\to \infty} \sup_{\theta\in \theta}\left|\frac{\partial F_{nj}(\theta,w)}{\partial \theta_i} - J_{i,j}(\theta)\right| = 0
\end{equation*}
for all $1\leq i\leq s$ and $1\leq j\leq s$, or, in other words, there exists a set $\mathcal{A}_1$ of probability $1$ such that
\begin{equation*}\lim_{n\to \infty} \frac{\partial F_n(\theta,w)}{\partial \theta_i} = J(\theta)\mbox{ uniformly on } \theta\mbox{ for all }w\in \mathcal{A}_1
\end{equation*}
Combining the obtained results with the fact that $J(\theta_0)$ is invertible by the hypothesis $(B)$ it follows that, for each fixed $w\in \mathcal{A}_1\cap \mathcal{A}_2$, letting $G(\theta)=J(\theta)$ for $\theta\in \Theta_0$, the sequence of functions $F_n(\theta,w)$ satisfies Theorem 2.1 and thus there exists $N>0$ and a sequence $(\theta_n(w))_{n\geq N}\in \theta$ such that \begin{equation*}
F_n(\theta_n(w),w)=0\mbox{ and }\lim_{n\to \infty} \theta_n(w) = \theta_0\mbox{ for all }w\in \mathcal{A}_1\cap \mathcal{A}_2
\end{equation*}
.

Thus $\hat{\theta}_n(X(w))=\theta_n(w)$ for all $w\in \mathcal{A}_1\cap \mathcal{A}_2\cap\mathcal{A}_3$ and it follows that 
\begin{equation*}
\lim_{n\to \infty} \hat{\theta}_n(X(w))=\lim_{n\to \infty} \theta_n(w)=\theta_0\mbox{ for all }w\in \mathcal{A}_1\cap \mathcal{A}_2\cap \mathcal{A}_3
\end{equation*}
Since the set $\mathcal{A}_1\cap \mathcal{A}_2\cap \mathcal{A}_3$ is a set of probability $1$ in $\Omega$, this proves that $\hat{\theta}_n(X)$ is a consistent estimator for $\theta_0$.

\vspace{0.3cm}

Asymptotic normality:
\vspace{0.3cm}

Once again denoting $(\beta_1,\cdots,\beta_s)=(\theta_1,\cdots,\theta_{s-r},\alpha_1,\cdots,\alpha_r)$, for each $i$ with $1\leq i\leq s$, and $w\in \mathcal{A}_1\cap \mathcal{A}_2\cap \mathcal{A}_3$, from the mean value theorem there must exists $t_i\in [0,1]$ and $y_i=t_i\theta_0+(1-t_i)\theta_n(w)$ such that
\begin{equation*}
\begin{aligned}
\frac{\partial}{\partial \beta_i} l\left(\theta_n(w),\alpha_0\ ;\ t\right) = \frac{\partial}{\partial \beta_i} l\left(\theta_0,\alpha_0\ ;\ t\right) + \sum_{j=1}^s  \frac{\partial l(y_i,\alpha_0\ ;\ 
t)}{\partial \theta_j \partial \beta_i}(\theta_{nj}(w) - \theta_{0j})\Rightarrow\\
\sum_{j=1}^s \frac{\partial l(y_i,\alpha_0\ ;\ 
t)}{\partial \theta_j \partial \beta_i}(\theta_{nj}(w) - \theta_{0j})=-\frac{\partial}{\partial \beta_i} l\left(\theta_0,\alpha_0\ ;\ t\right)
\end{aligned}
\end{equation*}
Thus, letting $A_n\in M_s(\R)$ be defined by $A_n=(a_{ij})$ where $a_{ij}=\frac{\partial l(y_i,\alpha_0\ ;\ 
t)}{\partial \theta_k \partial \beta_i}$ for $1\leq i\leq s$ and $1\leq i\leq s$, the above equation implies in
\begin{equation*} A_n (\theta_n(w)-\theta_0) = - v_n
\end{equation*}
where $v_n= \left(\frac{\partial}{\partial \beta_1} l\left(\theta^0,\alpha^0\ ;\ t\right),\cdots,\frac{\partial}{\partial \beta_s} l\left(\theta^0,\alpha^0\ ;\ t\right)\right)$

Now notice that $A_n$ converges in probability to $J(\theta_0)$, and from the Central Limit Theorem $v_n$ converges in distribution to $N(0,J(\theta_0))$. Thus, from the theorem of Lehman it follows that $\Theta_n(w)-\theta_0$ converges in distribution to
\begin{equation*} J(\theta_0)^{-1}N(0,J(\theta_0))=N(0,J(\theta_0)^{-1}),
\end{equation*}
which concludes the proof.
\end{proof}

\end{appendix}


\bibliographystyle{chicago}

\bibliography{reference}

\end{document}

