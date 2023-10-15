#include <Rcpp.h>

using namespace Rcpp;

// indexing(): Convert 2d index to 1d index to extract elements from dist objects

int indexing(int N, int i, int j){

  if(i == j){
    return NumericVector::get_na();
  }

  if (i < j+1){
    int temp;
    temp = i;
    i = j;
    j = temp;
  }

  return (2*N-1-j)*j/2 - 1 +(i-j);
}

//BSCalc(): Compute the two largest numbers in an array and their indexes given in ldiC (k=3).+
void BSCalc(NumericVector diC, IntegerVector ldiC, int NdiC, double &b, double &s, int &lb, int &ls){

  //temporary value
  double t;
  //temporary label
  int tl;

  b = diC[0];
  s = diC[1];

  lb = ldiC[0];
  ls = ldiC[1];

  //space for readability
  if(b > s){
    //switch(b,s)
    t = b;
    b = s;
    s = t;
    //switch(lb, ls)
    tl = lb;
    lb = ls;
    ls = tl;
  }

  for(int i = 2; i < NdiC; i++){

    //space for readability
    if(diC[i] < b){
      //modify b,s
      t = b;
      b = diC[i];
      s = t;
      //modify lb,ls
      tl = lb;
      lb = ldiC[i];
      ls = tl;
    } else if(diC[i] < s){ //modify s
      s = diC[i];
      //modify ls
      ls = ldiC[i];
    }
  }
}

//BSHCalc(): Compute the three largest numbers in an array and their indexes given in ldiC (k>=3).
void BSHCalc(NumericVector diC, IntegerVector ldiC, int NdiC, double &b, double &s, double &h, int &lb, int &ls, int &lh){

  //temporary values
  double t;
  double t2;
  //temporary labels
  int tl;
  int tl2;

  b = diC[0];
  s = diC[1];
  h = diC[2];

  lb = ldiC[0];
  ls = ldiC[1];
  lh = ldiC[2];

  if(b > s){//swap(b,s)
    t = b;
    b = s;
    s = t;
    //swap(ls,lb)
    tl = lb;
    lb = ls;
    ls = tl;
  }

  if(b > h){//swap(b,h)
    t = b;
    b = h;
    h = t;
    //swap(lb,lh)
    tl = lb;
    lb = lh;
    lh = tl;
  }

  if(s > h){//swap(s,h)
    t = s;
    s = h;
    h = t;
    //swap(ls,lh)
    tl = lh;
    ls = lh;
    ls = tl;
  }

  for(int i = 3; i < NdiC;i++){

    if(diC[i] < b){
      //modify b,s,h
      t = b;
      b = diC[i];
      t2 = s;
      s = t;
      h = t2;
      //modify lb,ls,lh
      tl = lb;
      lb = ldiC[i];
      tl2 = ls;
      ls = tl;
      lh = tl2;
    } else if(diC[i] < s){
      //modify s, h
      t = s;
      s = diC[i];
      h = t;
      //modify ls, lh
      tl = ls;
      ls = ldiC[i];
      lh = tl;
    } else if(diC[i] > h){
      //modify h
      h = diC[i];
      //modify lh
      lh = ldiC[i];
    }
  }
}

// UCT(): Compute Unit-Cluster Total distance matrices
NumericVector UCT(IntegerVector  C, NumericVector dist, int N, int k){

  NumericVector UCTdist(N*k);
  double dij;

  for(int i = 0; i < N; i++){
    for(int j = (i+1); j < N; j++){
      dij = dist[indexing(N, i,j)];
      UCTdist[i*k+C[j]] += dij;
      UCTdist[j*k+C[i]] += dij;
    }
  }
  return UCTdist;
}

//UCT2ASW: Compute the ASW from a UCT distance matrix
double UCT2ASW(IntegerVector C, NumericVector UCTdist, int N, int k, IntegerVector Nj, NumericVector ai, NumericVector bi, NumericVector si, NumericVector hi,
               IntegerVector lbi, IntegerVector lsi, IntegerVector lhi){

  double ASW = 0.;
  int li;
  int ik;
  NumericVector diC(k-1);
  IntegerVector ldiC(k-1);

  for(int i = 0; i < N; i++){

    ik = i*k;
    li = C[i];
    int l = 0;
    for(int j = 0; j < k; j++){
      if(j != li){
        diC[l] = UCTdist[ik+j]/Nj[j];
        ldiC[l] = j;
        l++;
      }
    }
    //space for readability

    if(k == 2){
      bi[i] = diC[0];
    } else if(k == 3){
      BSCalc(diC, ldiC, k-1, bi[i], si[i], lbi[i], lsi[i]);
    } else{
      BSHCalc(diC, ldiC, k-1, bi[i], si[i], hi[i], lbi[i], lsi[i], lhi[i]);
    }

    if(Nj[li] == 1){
      ai[i] = NumericVector::get_na();
      continue;
    } else{
      ai[i] = UCTdist[ik+li]/(Nj[li] - 1);
    }
    ASW += (bi[i] - ai[i])/std::max(bi[i], ai[i]);
  }
  return ASW/N;
}

//UpdateASW(): Efficiently compute the ASW in each swap
double UpdateASW(NumericVector UCTdist, IntegerVector C, NumericVector dist, int N, int i, int j, int k, IntegerVector Nj, NumericVector ai, NumericVector bi, NumericVector si, NumericVector hi,
                 IntegerVector lbi, IntegerVector lsi, IntegerVector lhi){
  IntegerVector updC = clone(C);
  updC[i] = j;
  double a;
  double b;
  double ASW = 0.;
  double td;
  int lk;
  int l2;
  double temp;  //temporary value
  double temp2;  //temporary value


  if(k == 2){ //Updating formulas for k == 2

    int l1 = abs(1-j);
    for(int l = 0; l < N; l++){
      lk = l*k;
      int l2 = updC[l];
      td = dist[indexing(N, l, i)];
      if(i == l){
        a = UCTdist[lk+j]/Nj[j];
        b = UCTdist[lk+l1]/(Nj[l1]-1);
        ASW += (b - a)/std::max(b, a);
      } else if(l2 == j){
        a = (UCTdist[lk+j] + td)/Nj[j];
        b = (UCTdist[lk+l1] - td)/(Nj[l1]-1);
        ASW += (b - a)/std::max(b, a);
      }  else{
        if( Nj[l1] == 2){
          ASW += 0;
        } else{
          a = (UCTdist[lk+l1] - td)/(Nj[l1] - 2);
          b = (UCTdist[lk+j] + td)/(Nj[j]+1);
          ASW += (b - a)/std::max(b, a);
        }
      }
    }
  } else if(k > 2){

    //space for readability
    int l1 = C[i];
    for(int l = 0; l < N; l++){

      lk = l*k;
      l2 = updC[l];
      td = dist[indexing(N, l, i)];
      if(l == i){ //x_l in Cl1; x_l = x_i
        a = UCTdist[lk+j]/Nj[j];

        if(j == lbi[l]){
          b = std::min(ai[l], si[l]);
        } else{
          b = std::min(ai[l], bi[l]);
        }
        ASW += (b - a)/std::max(b, a);
        //space for readability
      } else if(l2 == j){ //x_l in Cj; x_l != x_i

        a = (UCTdist[lk+j] + td)/Nj[j];
        temp = (UCTdist[lk+l1] - td)/(Nj[l1] - 1);

        if(l1 == lbi[l]){
          b = std::min(temp, si[l]);
        } else{
          b = std::min(temp, bi[l]);
        }
        ASW += (b - a)/std::max(b, a);
        //space for readability
      } else if(l2 != l1 && l2 != j){ //x_l not in Cj and Cl1
        if(Nj[l2] == 1){
          continue;
        } else{

          a = ai[l];
          temp = (UCTdist[lk+l1] - td)/(Nj[l1] - 1);
          temp2 = (UCTdist[lk+j] + td)/(Nj[j] + 1);

          if(k == 3){
            b = std::min(temp, temp2);
          } else{ //k > 3

            if(l1 != lbi[l] && j != lbi[l]){
              b = std::min({temp, temp2, bi[l]});
            } else{
              if(l1 == lsi[l] || j == lsi[l]){
                b = std::min({temp, temp2, hi[l]});
              } else{
                b = std::min({temp, temp2, si[l]});
              }
            }
          }

          ASW += (b - a)/std::max(b, a);
        }
        //space for readability
      } else { //x_l in Cl1 but x_l != x_i
        if( Nj[l1] == 2){
          continue;
        } else{
          a = (UCTdist[lk+l1] - td)/(Nj[l1] - 2);
          temp = (UCTdist[lk + j] + td)/(Nj[j] + 1);
          if(j == lbi[l]){
            b = std::min(temp, si[l]);
          } else{
            b = std::min(temp, bi[l]);
          }
          ASW += (b - a)/std::max(b, a);
        }
      }
    }

  }

  return ASW/N;
}

//UpdateUCT(): Update the UCT after finding the best swap
void UpdateUCT(NumericVector UCTdist, NumericVector dist, int N, int k, int i, int j, int li){

  double td;
  int lk;
  for(int l = 0; l < N; l++){
    lk = l*k;
    if(i == l) continue;
    td = dist[indexing(N, i,l)];
    UCTdist[lk + li] -= td;
    UCTdist[lk + j] += td;
  }
}

//UpdateUCT(): Update a,b,s terms after finding the best swap (k=3)
void UpdateABS(IntegerVector C, NumericVector UCTdist, int N, int k, IntegerVector Nj, NumericVector ai, NumericVector bi, NumericVector si,
               IntegerVector lbi, IntegerVector lsi){

  int li;
  int ik;
  NumericVector diC(k-1);
  IntegerVector ldiC(k-1);

  for(int i = 0; i < N; i++){

    ik = i*k;
    li = C[i];
    int l = 0;
    for(int j = 0; j < k; j++){
      if(j != li){
        diC[l] = UCTdist[ik+j]/Nj[j];
        ldiC[l] = j;
        l++;
      }
    }
    //space for readability

    BSCalc(diC, ldiC, k-1, bi[i], si[i], lbi[i], lsi[i]);


    if(Nj[li] == 1){
      ai[i] = NumericVector::get_na();
      continue;
    } else{
      ai[i] = UCTdist[ik+li]/(Nj[li] - 1);
    }
  }
}

//UpdateUCT(): Update a,b,s,h terms after finding the best swap (k>3)
void UpdateABSH(IntegerVector C, NumericVector UCTdist, int N, int k, IntegerVector Nj, NumericVector ai, NumericVector bi, NumericVector si, NumericVector hi,
                IntegerVector lbi, IntegerVector lsi, IntegerVector lhi){

  int li;
  int ik;
  NumericVector diC(k-1);
  IntegerVector ldiC(k-1);

  for(int i = 0; i < N; i++){

    ik = i*k;
    li = C[i];
    int l = 0;
    for(int j = 0; j < k; j++){
      if(j != li){
        diC[l] = UCTdist[ik+j]/Nj[j];
        ldiC[l] = j;
        l++;
      }
    }
    //space for readability

    BSHCalc(diC, ldiC, k-1, bi[i], si[i], hi[i], lbi[i], lsi[i], lhi[i]);

    if(Nj[li] == 1){
      ai[i] = NumericVector::get_na();
      continue;
    } else{
      ai[i] = UCTdist[ik+li]/(Nj[li] - 1);
    }
  }
}


//ASWOSil(): Compute the ASW (used inside OSil)
double ASWOSil(IntegerVector  C, NumericVector dist, int N, int k, IntegerVector Nj){

  NumericVector UCTdist(N*k);
  int li;
  double dij;

  double ASW_value = 0.;
  double a;
  double b;
  NumericVector diC(k-1);
  int ik;

  for(int i = 0; i < N; i++){
    for(int j = (i+1); j < N; j++){
      dij = dist[indexing(N, i,j)];
      UCTdist[i*k+C[j]] += dij;
      UCTdist[j*k+C[i]] += dij;
    }
  }

  for(int i = 0; i < N; i++){
    ik = i*k;
    li = C[i];

    if(Nj[li] == 1) continue;
    a = UCTdist[ik+li]/(Nj[li] - 1);

    int l = 0;
    for(int j = 0; j < k; j++){
      if(j != li){
        diC[l] = UCTdist[ik+j]/Nj[j];
        l++;
      }
    }

    b = min(diC);

    ASW_value += (b - a)/std::max(b, a);
  }
  return ASW_value/N;

}


//ASWCpp(): Compute the ASW
//' @export
// [[Rcpp::export(.ASWCpp)]]
double ASWCpp(IntegerVector C, NumericVector dist, int N, int k){

  IntegerVector Nj(k);
  NumericVector UCTdist(N*k);
  int li;
  double dij;

  double ASW_value = 0.;
  double a;
  double b;
  NumericVector diC(k-1);
  int ik;

  for(int i = 0; i < N; i++){ //efficient
    for(int j = 0; j < k; j++){
      if(C[i] == j){
        Nj[j] ++;
        break;
      }
    }
  }

  for(int i = 0; i < N; i++){
    for(int j = (i+1); j < N; j++){
      dij = dist[indexing(N, i,j)];
      UCTdist[i*k+C[j]] += dij;
      UCTdist[j*k+C[i]] += dij;
    }
  }

  for(int i = 0; i < N; i++){
    ik = i*k;
    li = C[i];

    if(Nj[li] == 1) continue;
    a = UCTdist[ik+li]/(Nj[li] - 1);

    int l = 0;
    for(int j = 0; j < k; j++){
      if(j != li){
        diC[l] = UCTdist[ik+j]/Nj[j];
        l++;
      }
    }

    b = min(diC);

    ASW_value += (b - a)/std::max(b, a);
  }
  return ASW_value/N;

}

//effOSilCpp(): Efficient Optimum Silhouette Clustering Algorithm
//' @export
// [[Rcpp::export(.effOSilCpp)]]
List effOSilCpp(NumericVector dist, IntegerVector iC, int N, int k){

  IntegerVector initC = clone(iC);
  IntegerVector Nj(k);
  IntegerVector Njtemp(k);
  NumericVector UCTdist(N*k);
  double bestASW;
  double tempASW;
  int li;

  NumericVector ai(N);
  NumericVector bi(N);
  NumericVector si(N);
  NumericVector hi(N);
  IntegerVector lbi(N);
  IntegerVector lsi(N);
  IntegerVector lhi(N);


  for(int i = 0; i < N; i++){ //efficient
    for(int j = 0; j < k; j++){
      if(initC[i] == j){
        Nj[j] ++;
        break;
      }
    }
  }

  int iter = 0;
  UCTdist = UCT(initC, dist, N, k);
  bestASW = UCT2ASW(initC, UCTdist, N, k, Nj, ai, bi, si, hi, lbi, lsi, lhi);

  bool swap;
  do{
    iter++;

    swap = false;
    IntegerVector swappedPair(2);

    for(int i = 0; i < N; i++){
      li = initC[i];

      if(Nj[li] == 1){
        continue;
      }

      for(int j = 0; j < k; j++){

        if(li == j){
          continue;
        }

        IntegerVector tempC = clone(initC);
        tempC[i] = j;
        tempASW = UpdateASW(UCTdist, initC, dist, N, i, j, k, Nj, ai, bi, si, hi, lbi, lsi, lhi);

        if(tempASW > bestASW){
          bestASW = tempASW;
          swappedPair[0] = i;
          swappedPair[1] = j;
          swap = true;
        }
      }
    }

    if(swap){
      li = initC[swappedPair[0]];
      Nj[li] --;
      Nj[swappedPair[1]] ++;
      initC[swappedPair[0]] = swappedPair[1];
      UpdateUCT(UCTdist, dist, N, k,
                swappedPair[0], swappedPair[1], li);
      if(k == 3){
        UpdateABS(initC, UCTdist, N, k, Nj, ai, bi, si, lbi, lsi);
      } else if(k > 3){
        UpdateABSH(initC, UCTdist, N, k, Nj, ai, bi, si, hi, lbi, lsi, lhi);
      }
    }
  }
  while (swap);
  return List::create(Named("Clustering") = initC+1L , _["ASW"] = bestASW,  _["nIter"] = iter);
}

//OSilCpp(): The original Optimum Silhouette Clustering Algorithm
//' @export
// [[Rcpp::export(.OSilCpp)]]
List OSilCpp(NumericVector dist, IntegerVector iC, int N, int k){

  IntegerVector initC = clone(iC);
  IntegerVector tempC(N);
  IntegerVector Nj(k);
  IntegerVector Njtemp(k);
  double bestASW;
  double tempASW;
  NumericVector UCTdist;
  int li;

  for(int i = 0; i < N; i++){ //efficient
    for(int j = 0; j < k; j++){
      if(initC[i] == j){
        Nj[j] ++;
        break;
      }
    }
  }

  int iter = 0;
  bestASW = ASWOSil(initC, dist, N, k, Nj);
  bool swap;
  do{
    iter++;

    swap = false;
    IntegerVector swappedPair(2);

    for(int i = 0; i < N; i++){
      li = initC[i];
      if(Nj[li] == 1){
        continue;
      }

      for(int j = 0; j < k; j++){

        if(li == j){
          continue;
        }

        Njtemp = clone(Nj);
        Njtemp[li] --;
        Njtemp[j] ++;

        IntegerVector tempC = clone(initC);
        tempC[i] = j;

        tempASW = ASWOSil(tempC, dist, N, k, Njtemp);

        if(tempASW > bestASW){

          bestASW = tempASW;
          swappedPair[0] = i;
          swappedPair[1] = j;
          swap = true;
        }
      }
    }

    if(swap){

      li = initC[swappedPair[0]];
      Nj[li] --;
      Nj[swappedPair[1]] ++;
      initC[swappedPair[0]] = swappedPair[1];
    }
  }
  while (swap);

  return List::create(Named("Clustering") = initC+1L, _["ASW"] = bestASW,  _["nIter"] = iter);
}
