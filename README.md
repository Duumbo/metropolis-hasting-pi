# Petit calcul de pi Monte-Carlo
Ce petit code calcule pi par chaine de Markov. Je t'envoie ce code pour montrer
l'importance du nombre d'états accessibles lors de la proposition d'un état dans
la chaine de Markov. La section 3.2 de mon mémoire (Chaine de Markov) se voulait
une brève explication de ce problème exactement, où l'on fait une erreur sur
la distribution de probabilité si l'on ne considère pas le nombre d'états accessibles
à un certain $x_\mu$. J'ai bien vu que ma section n'était pas claire, alors j'ai
réécrit certains passages pour mettre l'emphase là-dessus. C'est important car
nous n'avons pas pris ça en compte dans le code, et mVMC non plus. Dans la limite
où le nombre de site est très grand, c'est moins important car en moyenne chaque
état verra le même nombre d'état accessible, la probabilité de proposition d'un
état sera donc à peut-près uniforme, mais dans l'autre limite (ou même juste à
demi-rempli, où certains états bloquent d'autres par l'exclusion de Pauli), je
crois que ce code montre que c'est tout de même important. C'est d'autant plus
important pour les conditions de frontière fermées.

C'est aussi important de considérer ça lorsqu'on propose un spin flip, où que
l'on ne le propose pas avec la même probabilité que les autres.

# Comment ça marche
Une grille carrée de dimension SIZE, avec un cercle inscrit de rayon SIZE/2. On
génère une chaine de Markov pour marcher à travers la grille pour mesure le
ratio entre l'aire du cercle et du carrée. Le seul point important si tu veux jouer
avec c'est la variable NEIGHBORS_CHECK: à false les états sont acceptés toujours
(distribution de probabilité uniforme avec le bilan détaillé)
à true on y multiplie le ratio du nombre d'état accessible entre le
précédent et le nouveau.

Quelques stats que j'ai obtenu avec ça
100 runs de
const SIZE: usize = 100;
const NSAMP: usize = 10_000_000;
const NWARM: usize = 10_000;

Neighbors check
Avg = 3.1392757999999996
Dev = 0.00048130196640627787

No Neighbors check
Avg = 3.1810236959999996
Dev = 0.0005005803373911988

# Run
```bash
cargo run -r --bin mc-prop-acc
```

# Snippet mVMC relevant
Dans cette portion de code, ils devraient prendre en considération la même chose
que je parle ici, mais ils ne le font pas. Ça devrait être
`p = w * (n_accessible / n_previous_accessible)`, pour que ça fonctionne avec le
cas à conditions de frontière fermées.
```C
void VMCMakeSample(...)
{
    (...)
          if(updateType==HOPPING) { /* hopping */
        Counter[0]++;

        StartTimer(31);
        makeCandidate_hopping(&mi, &ri, &rj, &s, &rejectFlag,
                              TmpEleIdx, TmpEleCfg);
        StopTimer(31);

        if(rejectFlag) continue;

        StartTimer(32);
        StartTimer(60);
        /* The mi-th electron with spin s hops to site rj */
        updateEleConfig(mi,ri,rj,s,TmpEleIdx,TmpEleCfg,TmpEleNum);
        UpdateProjCnt(ri,rj,s,projCntNew,TmpEleProjCnt,TmpEleNum);
    (...)
            /* Metroplis */
        x = LogProjRatio(projCntNew,TmpEleProjCnt);
        w = exp(2.0*(x+creal(logIpNew-logIpOld)));
        if( !isfinite(w) ) w = -1.0; /* should be rejected */

        if(w > genrand_real2()) { /* accept */
          // UpdateMAll will change SlaterElm, InvM (including PfM)
          StartTimer(63);
          UpdateMAll(mi,s,TmpEleIdx,qpStart,qpEnd);
          //            UpdateMAll(mi,s,TmpEleIdx,qpStart,qpEnd);
          StopTimer(63);

          for(i=0;i<NProj;i++) TmpEleProjCnt[i] = projCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          Counter[1]++;
        } else { /* reject */
          revertEleConfig(mi,ri,rj,s,TmpEleIdx,TmpEleCfg,TmpEleNum);
        }
        StopTimer(32);

      } else if(updateType==EXCHANGE) { /* exchange */

}

void makeCandidate_hopping(int *mi_, int *ri_, int *rj_, int *s_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg) {
  const int icnt_max = Nsite*Nsite;
  int icnt;
  int mi, ri, rj, s, flag;

  flag = 0; // FALSE
  do {
    mi = gen_rand32()%Ne;
    s = (genrand_real2()<0.5) ? 0 : 1;
    ri = eleIdx[mi+s*Ne];
  } while (LocSpn[ri] == 1);

  icnt = 0;
  do {
    rj = gen_rand32()%Nsite;
    if(icnt> icnt_max){
      flag = 1; // TRUE
      break;
    }
    icnt+=1;
  } while (eleCfg[rj+s*Nsite] != -1 || LocSpn[rj]==1);

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_ = s;
  *rejectFlag_ = flag;

  return;
}
```
