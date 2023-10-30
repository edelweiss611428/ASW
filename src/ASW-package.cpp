#include <Rcpp.h>

using namespace Rcpp;

// 0. Relevant functions
// indexing(): Convert 2d index to 1d index to extract elements from dist objects

int indexing(int N, int i, int j){

  if(i == j){
    return NumericVector::get_na();
  }

  if (i < j){
    int temp;
    temp = i;
    i = j;
    j = temp;
  }

  return ((2*N-1-j)*j >> 1) - 1 +(i-j);
}


//ASWCpp(): Compute the Average Silhouette Width

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

  int iij;

  for(int i = 0; i < N; i++){
    iij = indexing(N, i,i+1);
    for(int j = (i+1); j < N; j++){
      dij = dist[iij];
      UCTdist[i*k+C[j]] += dij;
      UCTdist[j*k+C[i]] += dij;
      iij++;
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


// subDistCpp(): Extract a sub-distance matrix from a "dist" object

// [[Rcpp::export(.subDistCpp)]]
NumericVector subDistCpp(NumericVector dist, IntegerVector idx, bool diag, bool upper, int N, int n){

  NumericVector subDMat((n-1)*n >> 1);
  int l = 0;

  for(int i = 0; i < n-1; i++){
    for(int j = (i+1); j < n; j++){
      if(idx[i] == idx[j]){
        subDMat[l] = 0;
      } else{
        subDMat[l] = dist[indexing(N, idx[i], idx[j])];
      }
      l++;
    }
  }

  subDMat.attr("Size") = n;
  subDMat.attr("Diag") = diag;
  subDMat.attr("Upper") = upper;
  subDMat.attr("class") = "dist";

  return subDMat;
}


//SWCpp(): Compute the Silhouette Widths

// [[Rcpp::export(.SWCpp)]]
List SWCpp(IntegerVector C, NumericVector dist, int N, int k){

  IntegerVector Nj(k);
  NumericVector UCTdist(N*k);
  int li;
  double dij;

  NumericVector sil_width(N);
  double a;
  double b;
  IntegerVector neighbor(N);
  NumericVector diC(k-1);
  IntegerVector ldiC(k-1);
  int ik;
  int l;

  for(int i = 0; i < N; i++){ //efficient
    for(int j = 0; j < k; j++){
      if(C[i] == j){
        Nj[j] ++;
        break;
      }
    }
  }

  int iij;

  for(int i = 0; i < N; i++){
    iij = indexing(N, i,i+1);
    for(int j = (i+1); j < N; j++){
      dij = dist[iij];
      UCTdist[i*k+C[j]] += dij;
      UCTdist[j*k+C[i]] += dij;
      iij++;
    }
  }

  for(int i = 0; i < N; i++){
    ik = i*k;
    li = C[i];

    if(Nj[li] == 1) continue;
    a = UCTdist[ik+li]/(Nj[li] - 1);

    l = 0;
    for(int j = 0; j < k; j++){
      if(j != li){
        diC[l] = UCTdist[ik+j]/Nj[j];
        ldiC[l]=j;
        l++;
      }
    }

    l = which_min(diC);
    b = diC[l];
    neighbor[i] = ldiC[l];
    sil_width[i] = (b - a)/std::max(b, a);
  }

  return List::create(Named("neighbor") = neighbor, _["sil_width"] = sil_width);

}



// 1. The Efficient Optimum Silhouette algorithm

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
  int iij;
  for(int i = 0; i < N; i++){
    iij = indexing(N, i,i+1);
    for(int j = (i+1); j < N; j++){
      dij = dist[iij];
      UCTdist[i*k+C[j]] += dij;
      UCTdist[j*k+C[i]] += dij;
      iij++;
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


//UpdateABSH(): Update a,b,s,h terms after finding the best swap (k>=3)
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
    if(k == 3){
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
  }
}


//effOSilCpp(): the Efficient Optimum Silhouette Clustering Algorithm

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
      if(k >= 3){
        UpdateABSH(initC, UCTdist, N, k, Nj, ai, bi, si, hi, lbi, lsi, lhi);
      }
    }
  }
  while (swap);
  return List::create(Named("Clustering") = initC+1L , _["ASW"] = bestASW,  _["nIter"] = iter);
}


// 2. The original Optimum Silhouette algorithm

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
  int iij;

  for(int i = 0; i < N; i++){
    iij = indexing(N, i,i+1);
    for(int j = (i+1); j < N; j++){
      dij = dist[iij];
      UCTdist[i*k+C[j]] += dij;
      UCTdist[j*k+C[i]] += dij;
      iij++;
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


//OSilCpp(): The original Optimum Silhouette Clustering Algorithm
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




// The Scalable Optimum Silhouette algorithm

//scalOSil_PC_Step() performs the PC step of scalOSil

// [[Rcpp::export(.scalOSil_PC)]]
List scalOSil_PC_Step(NumericVector dist, IntegerVector iC, int N, int k){

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
      if(k >= 3){
        UpdateABSH(initC, UCTdist, N, k, Nj, ai, bi, si, hi, lbi, lsi, lhi);
      }
    }
  }
  while (swap);

  int ik;
  int l2;

  if(k == 2){

    for(int i = 0; i < N; i++){
      ik = i << 1;
      li = initC[i];
      l2 = abs(1-li);

      if(Nj[li] == 1){
        ai[i] = NumericVector::get_na();
        bi[i] = UCTdist[ik+l2]/Nj[l2];
      } else{
        ai[i] = UCTdist[ik+li]/(Nj[li]-1);
        bi[i] = UCTdist[ik+l2]/Nj[l2];
      }

    }

  }

  return List::create(Named("PC") = initC, _["ASW"] = bestASW,  _["PC_UCT"] = UCTdist, _["Nj"] = Nj,
                      _["PC_ai"] = ai, _["PC_bi"] = bi, _["PC_si"] = si, _["PC_lbi"] = lbi);
}

//miniSW() computes the SW only for the newest unassigned point
double miniSW(NumericVector diC, int m, int k){

  NumericVector TempdiC(k-1);
  int l = 0;
  for(int i = 0; i < k; i++){
    if(i!=m){
      TempdiC[l] = diC[i];
      l++;
    }
  }
  double a = diC[m];
  double b = min(TempdiC);
  return (b-a)/std::max(a,b);

}

//scalOSil_C_Step() performs the C step of scalOSil

// [[Rcpp::export(.scalOSil_C)]]
IntegerVector scalOSil_C_Step(NumericVector dist, int k, List PC_result, IntegerVector idxPC, IntegerVector idxC, int n1, int n2, int N){

  //Copy values from PC_result

  IntegerVector PC = PC_result["PC"];
  NumericVector PC_UCT = PC_result["PC_UCT"];
  IntegerVector Nj = PC_result["Nj"];
  NumericVector PC_ai = PC_result["PC_ai"];
  NumericVector PC_bi = PC_result["PC_bi"];
  NumericVector PC_si = PC_result["PC_si"];
  IntegerVector PC_lbi = PC_result["PC_lbi"];


  //Final clustering

  IntegerVector FC(N);

  //Copy values from PC to FC;

  for(int i = 0; i < n1; i++){
    FC[i] = PC[i];
  }

  //Set initial value for counter l

  int l = n1;

  //Index for each unassigned unit

  int ii;

  //Others

  int lj;
  int ij; //index for j when updating diC
  double a;
  double b;
  double td;
  double tempASW;
  double bestASW;

  double t;

  NumericVector TD(n1); //Needed for efficiency
  NumericVector diC(k);


  for(int i = 0; i < n2; i++){ //Unassigned units

    NumericVector TdiC(k);

    ii = idxC[i];

    //Update diC

    for(int j = 0; j < n1; j++){

      lj = PC[j];
      ij = idxPC[j];

      td = dist[indexing(N, ii, ij)];
      TD[j] = td;
      TdiC[lj] += td;

    }

    for(int j = 0; j < k; j++){
      diC[j] = TdiC[j]/Nj[j];
    }

    //Consider each assignment

    bestASW = -1;

    for(int m = 0; m < k; m++){

      tempASW = 0;

      //Efficiently compute "tempASW": First, for all units in PC, and then, for the unassigned unit i

      for(int j = 0; j < n1; j++){

        lj = PC[j];

        if(lj == m){

          if(Nj[lj] == 1){
            a = TD[j];
            b = PC_bi[j];
          } else{

            a = (PC_UCT[j*k+m] + TD[j])/Nj[lj];
            b = PC_bi[j];

          }

          tempASW += (b - a)/std::max(b, a);

        } else{

          if(Nj[lj] == 1){
            continue;
          } else{

            a = PC_ai[j];
            t = (PC_UCT[j*k+m] + TD[j])/(Nj[m] + 1);

            if(k == 2){

              b = t;

            } else{

              if(m != PC_lbi[j]){
                b = std::min(PC_bi[j], t);
              } else{
                b = std::min(PC_si[j], t);
              }

            }

            tempASW += (b - a)/std::max(b, a);
          }

        }

      } // end of for(int j = 0; j < n1; j++) //


      tempASW += miniSW(diC, m, k);
      tempASW = tempASW/(n1+1);

      if(tempASW > bestASW){

        bestASW = tempASW;
        FC[l] = m;

      }

    } //end of for(int m = 0; m < k; m++)//

    l++;

  }
  return FC;

}

// 4. The original Fast Optimum Silhouette Algorithm

//expand_dist() expands the dist object by one row
NumericVector expand_dist(NumericVector dist, int N){

  NumericVector expanded_dist((N+1)*N>>1);
  int l = 0;
  for(int j = 0; j < N ; j++){
    for(int i = (j+1); i < N+1; i++){
      if(i == N){
        l++;
        continue;
      }
      expanded_dist[l] = dist[indexing(N,i,j)];
      l++;
    }
  }

  expanded_dist.attr("Size") = N+1;
  expanded_dist.attr("Diag") = false;
  expanded_dist.attr("Upper") = false;
  expanded_dist.attr("class") = "dist";

  return expanded_dist;
}

//FOSil_PC_Step() performs the PC step of FOSil

//[[Rcpp::export(.FOSil_PC)]]
List FOSil_PC_Step(NumericVector dist, IntegerVector iC, int N, int k){

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

  return List::create(Named("PC") = initC, _["ASW"] = bestASW, _["Nj"] = Nj);
}

//FOSil_C_Step() performs the C step of FOSil

//[[Rcpp::export(.FOSil_C)]]
IntegerVector FOSil_C_Step(NumericVector dist, NumericVector distPC, int k, List PC_result, IntegerVector idxPC, IntegerVector idxC, int n1, int n2, int N){

  //Copy values from PC_result

  IntegerVector PC = PC_result["PC"];
  IntegerVector Nj = PC_result["Nj"];

  //Final clustering

  IntegerVector FC(N);

  //Others

  IntegerVector tempC(n1+1);
  IntegerVector Njtemp(k);
  double tempASW;
  double bestASW;
  int n3 = n1+1;
  int l = n1;

  //Copy values from PC to FC and tempC;

  for(int i = 0; i < n1; i++){
    FC[i] = PC[i];
    tempC[i] = PC[i];
  }



  NumericVector dist2 = expand_dist(distPC, n1);

  for(int i = 0; i < n2; i++){
    bestASW = -1;

    for(int j = 0; j < n1; j++){
      dist2[indexing(n3, n1,j)] = dist[indexing(N, idxPC[j], idxC[i])];
    }

    for(int j = 0; j < k; j++){
      tempC[n1] = j;
      Njtemp = clone(Nj);
      Njtemp[j]++;
      tempASW = ASWOSil(tempC, dist2, n1+1,k,Njtemp);
      if(tempASW > bestASW){
        bestASW = tempASW;
        FC[l] = j;
      }
    }
    l++;
  }


  return FC;

}

//5. The PAMSil algorithm


//nearest_and_2nd_nearest() computes the  nearest and second nearest distances from each point to other medoids.
void nearest_and_2nd_nearest(NumericVector dist, IntegerVector medoids, IntegerVector C1st, IntegerVector &C2nd, NumericVector &d1st, NumericVector &d2nd, int N, int k){

  NumericVector diC(k-1);
  IntegerVector ldiC(k-1);
  int l;

  int l1st;
  int idx_min;

  for(int i = 0; i < N; i++){
    l1st = C1st[i];
    l = 0;

    for(int j = 0; j < k; j++){

      if(l1st != j){

        diC[l] = dist[indexing(N, i, medoids[j])];
        ldiC[l] = j;
        l++;

      }

    }

    idx_min = which_min(diC);
    C2nd[i] = ldiC[idx_min];
    d2nd[i] = diC[idx_min];

    if(i == medoids[l1st]){
      d1st[i] = 0;
    } else {
      d1st[i] = dist[indexing(N, i, medoids[l1st])];
    }
  }
}

//swap_to_clustering() efficiently compute a clustering given a swap in O(N) time
IntegerVector swap_to_clustering(NumericVector dist, IntegerVector medoids, IntegerVector C1st, IntegerVector C2nd, NumericVector d1st, NumericVector d2nd,
                                 int N, int k, int s_i, int s_j){

  IntegerVector C_new(N);
  int li;
  double td;

  for(int i = 0; i < N; i++){
    li = C1st[i];

    if(i == s_i){
      C_new[i] = s_j;
      continue;
    } else if(i == medoids[li]){
      C_new[i] = li;
      continue;
    }

    td = dist[indexing(N,i,s_i)];

    if(li != s_j){

      if(td < d1st[i]){
        C_new[i] = s_j;
      } else{
        C_new[i] = li;
      }

    } else{

      if(td < d2nd[i]){
        C_new[i] = s_j;
      } else{
        C_new[i] = C2nd[i];
      }
    }

  }

  return C_new;

}

//PAMSilCpp() implements the PAMSil algorithm

// [[Rcpp::export(.PAMSilCpp)]]
List PAMSilCpp(NumericVector dist, IntegerVector C, IntegerVector medoids, int N, int k){

  IntegerVector fC = clone(C);
  IntegerVector fmedoids = clone(medoids);
  IntegerVector tempC(N);
  IntegerVector C2nd(N);
  NumericVector d1st(N);
  NumericVector d2nd(N);
  nearest_and_2nd_nearest(dist, fmedoids, fC, C2nd, d1st, d2nd, N, k);
  int iter = 0;
  bool swap;
  double tempASW;


  double bestASW = ASWCpp(fC, dist, N, k);

  do{
    iter++;

    swap = false;
    IntegerVector swappedPair(2);


    for(int i = 0; i < N; i++){
      if(i == fmedoids[fC[i]]) continue;
      for(int j = 0; j < k; j++){

        tempC = swap_to_clustering(dist, medoids, fC, C2nd, d1st, d2nd, N, k, i, j);
        tempASW = ASWCpp(tempC, dist, N, k);

        if(tempASW > bestASW){

          bestASW = tempASW;
          swappedPair[0] = j;
          swappedPair[1] = i;
          swap = true;

        }

      }

    }

    if(swap){

      fmedoids[swappedPair[0]] = swappedPair[1];
      fC = swap_to_clustering(dist, fmedoids, fC, C2nd, d1st, d2nd, N, k, swappedPair[1], swappedPair[0]);
      nearest_and_2nd_nearest(dist, fmedoids, fC, C2nd, d1st, d2nd, N, k);

    }

  } while (swap);

  return List::create(Named("Clustering") = fC+1, _["medoids"] = fmedoids+1, _["ASW"] = bestASW,  _["nIter"] = iter);

}

