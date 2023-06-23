#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define diziBoyut 20
#define emir 4
struct get_function
{
    double P1[diziBoyut];
    double P2[diziBoyut];
    double PTurev1[diziBoyut];
    double PTurev2[diziBoyut];
    double U1[diziBoyut];
    double U2[diziBoyut];
    double U3[diziBoyut];
    double U4[diziBoyut];
    double U5[diziBoyut];
    double L1[diziBoyut];
    double L2[diziBoyut];
    double L3[diziBoyut];
    double L4[diziBoyut];
    double L5[diziBoyut];
    double T1[diziBoyut];
    double T2[diziBoyut];
    double T3[diziBoyut];
    double T4[diziBoyut];
    double T5[diziBoyut];
    double TersT1[diziBoyut];
    double TersT2[diziBoyut];
    double TersT3[diziBoyut];
    double TersT4[diziBoyut];
    char fonkKontrol;
    int trigKontrol[diziBoyut];
    int tersTrigKontrol[diziBoyut];
    double Matris[diziBoyut][diziBoyut];
    int n;
};
    struct get_function function1;
void turevAl()
{
    int i=0;
    while(function1.P1[i]!='\0')
    {
      function1.PTurev1[i] = function1.P1[i] * function1.P2[i];
      function1.PTurev2[i] = function1.P2[i] - 1;
      if(function1.PTurev2[i]<0)
      {
          function1.PTurev2[i]=1;
      }
      i++;
    }
}
double turev_answer(double x)
{
    int i=0;
    double sonuc=0;
    while(function1.PTurev1[i] != '\0')
    {
        sonuc += function1.PTurev1[i]*(pow(x,function1.PTurev2[i]));
        i++;
    }
    return sonuc;
}
int Menu()
    {
        int choice;
        printf("Quit: 0\n");
        printf("Bisection: 1\n");
        printf("Regula-Falsi : 2\n");
        printf("Newton Raphson : 3\n");
        printf("Inverse Matrix : 4\n");
        printf("Gauss Elimination : 5\n");
        printf("Gauss-Seidel : 6\n");
        printf("Numerical Differentiation : 7\n");
        printf("Simpson's Rule : 8\n");
        printf("Trapezoidal Rule : 9\n");
        printf("Gregory-Newton : 10\n");

       printf("Choice: ");
       scanf("%d",&choice);
    return choice;
    }
double function_answer(double x)
{
    int i=0;
    double sum=0;
    while(function1.P1[i]!=0)
    {
        sum+=function1.P1[i]*(pow(x,function1.P2[i]));
        i++;
    }
    i=0;
    while(function1.U1[i]!=0)
    {
        sum+=(function1.U4[i]*pow(pow(function1.U5[i],function1.U1[i]*pow(x,function1.U2[i])),function1.U4[i]));
        i++;
    }
    i=0;
    while(function1.L1[i]!=0)
    {
        sum+=pow(log10(function1.L1[i]*pow(x,function1.L2[i]))/log10(function1.L5[i]),function1.L4[i])*function1.L3[i];
        i++;
    }
    i=0;
    while(function1.T1[i]!=0)
    {
            if(function1.trigKontrol[i]==0) // sin
            {
                sum += function1.T3[i] * (pow(sin(function1.T1[i] * pow(x,function1.T2[i])),function1.T4[i]));
            }
             else if(function1.trigKontrol[i]==1) //cos
            {
                sum += function1.T3[i] * (pow(cos(function1.T1[i] * pow(x,function1.T2[i])),function1.T4[i]));
            }
            else if(function1.trigKontrol[i]==2)  // tan
            {
                sum += function1.T3[i] * (pow(tan(function1.T1[i] * pow(x,function1.T2[i])),function1.T4[i]));
            }
            else if(function1.trigKontrol[i]==3)  //cot
            {
                 sum += function1.T3[i] * (pow(1/tan(function1.T1[i] * pow(x,function1.T2[i])),function1.T4[i]));
            }
        i++;
    }
    i=0;
    while(function1.TersT1[i]!=0)
    {
            if(function1.tersTrigKontrol[i]==0) // sin
            {
                sum += function1.TersT3[i] * (pow(asin(function1.TersT1[i] * pow(x,function1.TersT2[i])),function1.TersT4[i]));
            }
             else if(function1.tersTrigKontrol[i]==1) //cos
            {
                sum += function1.TersT3[i] * (pow(acos(function1.TersT1[i] * pow(x,function1.TersT2[i])),function1.TersT4[i]));
            }
            else if(function1.tersTrigKontrol[i]==2)  // tan
            {
                sum += function1.TersT3[i] * (pow(atan(function1.TersT1[i] * pow(x,function1.TersT2[i])),function1.TersT4[i]));
            }
            else if(function1.tersTrigKontrol[i]==3)  //cot
            {
                 sum += function1.TersT3[i] * (pow(1/atan(function1.TersT1[i] * pow(x,function1.TersT2[i])),function1.TersT4[i]));
            }
        i++;
    }
    return sum;
}
void fonk_al()
{
    int i=0;
    printf("polinom fonk. girmek icin a'ya ustel girmek icin e'ye basiniz: ");
    scanf(" %c",&function1.fonkKontrol);
    while(function1.fonkKontrol!='e'){ // polinom fonk. alýmý
    if(function1.fonkKontrol!='a' && function1.fonkKontrol!='e')
    {
        printf("lutfen gecerli bir tusa basiniz!!\n");
        printf("polinom fonk. girmek icin a'ya ustel girmek icin e'ye basiniz: ");
        scanf(" %c",&function1.fonkKontrol);
    }else
        {
    printf("x'in katsayisini girin: ");
    scanf("%lf",&function1.P1[i]);
    printf("x'in ussunu girin: ");
    scanf("%lf",&function1.P2[i]);
    printf("polinom fonk. girmek icin a'ya ustel girmek icin e'ye basiniz: ");
    scanf(" %c",&function1.fonkKontrol);
    printf("polinom katsayilar: \n%lf\n%lf\n\n",function1.P1[i],function1.P2[i]);
    i++;
        }
                             }
      i=0;
    if(function1.fonkKontrol=='e')
    {
    function1.fonkKontrol='a';
    printf("ustel fonk. girmek icin a'ya logaritmik girmek icin e'ye basiniz: ");
    scanf(" %c",&function1.fonkKontrol);
    while(function1.fonkKontrol!='e'){// ustel fonk. alýmý
    if(function1.fonkKontrol!='a' && function1.fonkKontrol!='e')
    {
        printf("lutfen gecerli bir tusa basiniz!!\n");
        printf("ustel fonk. girmek icin a'ya logaritmik girmek icin e'ye basiniz: ");
        scanf(" %c",&function1.fonkKontrol);
    }else
        {
    printf("x'in katsayisini girin: ");
    scanf("%lf",&function1.U1[i]);
    printf("x'in ussunu girin: ");
    scanf("%lf",&function1.U2[i]);
    printf("fonk. katsayi girin:  ");
    scanf("%lf",&function1.U3[i]);
    printf("fonk. ussunu girin:  ");
    scanf("%lf",&function1.U4[i]);
    printf("log tabani girin:  ");
    scanf("%lf",&function1.U5[i]);
    printf("ustel fonk. girmek icin a'ya logaritmik girmek icin e'ye basiniz: ");
    scanf(" %c",&function1.fonkKontrol);
    printf("ustel katsayilar: \n%lf\n%lf\n%lf\n%lf\n%lf\n\n",function1.U1[i],function1.U2[i],function1.U3[i],function1.U4[i],function1.U5[i]);
    i++;
        }
                             }
    }
    i=0;
    if(function1.fonkKontrol=='e')
    {
    function1.fonkKontrol='a';
    printf("logaritmik fonk. girmek icin a'ya trigonometrik girmek icin e'ye basiniz: ");
    scanf(" %c",&function1.fonkKontrol);
    while(function1.fonkKontrol!='e'){// logaritmik fonk alýmý
    if(function1.fonkKontrol!='a' && function1.fonkKontrol!='e')
    {
        printf("lutfen gecerli bir tusa basiniz!!\n");
        printf("logaritmik fonk. girmek icin a'ya trigonometrik girmek icin e'ye basiniz: ");
        scanf(" %c",&function1.fonkKontrol);
    }else
        {
    printf("x'in katsayisini girin: ");
    scanf("%lf",&function1.L1[i]);
    printf("x'in ussunu girin: ");
    scanf("%lf",&function1.L2[i]);
    printf("fonk. katsayi girin:  ");
    scanf("%lf",&function1.L3[i]);
    printf("fonk. ussunu girin:  ");
    scanf("%lf",&function1.L4[i]);
    printf("log tabani girin:  ");
    scanf("%lf",&function1.L5[i]);
    printf("logaritmik fonk. girmek icin a'ya trigonometrik girmek icin e'ye basiniz: ");
    scanf(" %c",&function1.fonkKontrol);
    printf("logaritmik katsayilar: \n%lf\n%lf\n%lf\n%lf\n%lf\n\n",function1.L1[i],function1.L2[i],function1.L3[i],function1.L4[i],function1.L5[i]);
    i++;
        }
                             }
    }
    i=0;
    if(function1.fonkKontrol=='e')
    {
    function1.fonkKontrol='a';
    printf("trigonometrik fonk. girmek icin a'ya ters trigonometrik girmek icin e'ye basiniz: ");
    scanf(" %c",&function1.fonkKontrol);
    while(function1.fonkKontrol!='e'){ // trigonometrik fonk alýmý
    if(function1.fonkKontrol!='a' && function1.fonkKontrol!='e')
    {
        printf("lutfen gecerli bir tusa basiniz!!\n");
        printf("trigonometrik fonk. girmek icin a'ya ters trigonometrik girmek icin e'ye basiniz: ");
        scanf(" %c",&function1.fonkKontrol);
    }else
        {
    printf(" sin: 0\n cos: 1\n tan:2\n cot: 3\n secim: ");
    scanf("%d",&function1.trigKontrol[i]);
    printf("x'in katsayisini girin: ");
    scanf("%lf",&function1.T1[i]);
    printf("x'in ussunu girin: ");
    scanf("%lf",&function1.T2[i]);
    printf("fonk. katsayi girin:  ");
    scanf("%lf",&function1.T3[i]);
    printf("fonk. ussunu girin:  ");
    scanf("%lf",&function1.T4[i]);
    printf("trigonometrik fonk. girmek icin a'ya ters trigonometrik girmek icin e'ye basiniz: ");
    scanf(" %c",&function1.fonkKontrol);
    printf("trigonometrik katsayilar: \n%lf\n%lf\n%lf\n%lf\n\n",function1.T1[i],function1.T2[i],function1.T3[i],function1.T4[i]);
    i++;
        }
                             }
    }
     i=0;
    if(function1.fonkKontrol=='e')
    {
    function1.fonkKontrol='a';
    printf("ters trigonometrik fonk. girmek icin a'ya cikmak girmek icin e'ye basiniz: ");
    scanf(" %c",&function1.fonkKontrol);
    while(function1.fonkKontrol!='e'){ // ters trigonometrik fonk alýmý
    if(function1.fonkKontrol!='a' && function1.fonkKontrol!='e')
    {
        printf("lutfen gecerli bir tusa basiniz!!\n");
        printf("ters trigonometrik fonk. girmek icin a'ya cikmak icin e'ye basiniz: ");
        scanf(" %c",&function1.fonkKontrol);
    }else
        {
    printf(" arcsin: 0\n arccos: 1\n arctan:2\n arccot: 3\n secim: ");
    scanf("%d",&function1.tersTrigKontrol[i]);
    printf("x'in katsayisini girin: ");
    scanf("%lf",&function1.TersT1[i]);
    printf("x'in ussunu girin: ");
    scanf("%lf",&function1.TersT2[i]);
    printf("fonk. katsayi girin:  ");
    scanf("%lf",&function1.TersT3[i]);
    printf("fonk. ussunu girin:  ");
    scanf("%lf",&function1.TersT4[i]);
    printf("ters trigonometrik fonk. girmek icin a'ya cikmak girmek icin e'ye basiniz: ");
    scanf(" %c",&function1.fonkKontrol);
    printf("ters trigonometrik katsayilar: \n%lf\n%lf\n%lf\n%lf\n\n",function1.TersT1[i],function1.TersT2[i],function1.TersT3[i],function1.TersT4[i]);
    i++;
        }
                             }
    }
}
void bisection(){
    double mid,start, end, tol;
    int i = 0,kontrol=1,max_iter,durmaKosulu;
    printf("baslangic degerini girin: ");
    scanf("%lf",&start);
    printf("son degerini girin: ");
    scanf("%lf",&end);
    printf("hata miktarini girin: ");
    scanf("%lf",&tol);
    printf("max iterasyonu girin: ");
    scanf("%d",&max_iter);
    printf("1. f(x) <= epsilon\n2. end-start/2 uzeri n <= epsilon\n");
    printf("Durma Kosulu Girin: ");
    scanf("%d",&durmaKosulu);
    if(durmaKosulu==1)
    {
            while (function_answer(mid) < tol && i < max_iter && kontrol==1)
        {
            i++;
            mid = (start + end) / 2;
            printf("start     : %lf\n",start);
            printf("end       : %lf\n",end);
            printf("mid       : %lf\n",mid);
            printf("f(start)  : %lf\n",function_answer(start));
            printf("f(end)    : %lf\n",function_answer(end));
            printf("f(mid)    : %lf\n",function_answer(mid));
            printf("iteration : %d\n",i);
            if (function_answer(mid) == 0) {
                kontrol=0;
            }
            else if (function_answer(mid) * function_answer(start) < 0) {
                end = mid;
            }
            else if(function_answer(mid) * function_answer(end) < 0)
            {
                start = mid;
            }
            else
                kontrol=0;
        }
    }
        i=0;
        if(durmaKosulu==2)
        {
            while ((end - start) / pow(2,i) >= tol && i < max_iter && kontrol==1)
        {
            i++;
            mid = (start + end) / 2;
            printf("start     : %lf\n",start);
            printf("end       : %lf\n",end);
            printf("mid       : %lf\n",mid);
            printf("f(start)  : %lf\n",function_answer(start));
            printf("f(end)    : %lf\n",function_answer(end));
            printf("f(mid)    : %lf\n",function_answer(mid));
            printf("iteration : %d\n",i);
            if (function_answer(mid) == 0) {
                kontrol=0;
            }
            else if (function_answer(mid) * function_answer(start) < 0) {
                end = mid;
            }
            else if(function_answer(mid) * function_answer(end) < 0)
            {
                start = mid;
            }
            else
                kontrol=0;
        }
    }
    printf("Kok: %lf\n", mid);
}
void regula_falsi()
{
    double start, end, tol,point;
    int i = 0,kontrol=1,max_iter,durmaKosulu;
    printf("baslangic degerini girin: ");
    scanf("%lf",&start);
    printf("son degerini girin: ");
    scanf("%lf",&end);
    printf("hata miktarini girin: ");
    scanf("%lf",&tol);
    printf("max iterasyonu girin: ");
    scanf("%d",&max_iter);
    printf("1. f(x) <= epsilon\n2. end-start/2 uzeri n <= epsilon\n");
    printf("Durma Kosulu Girin: ");
    scanf("%d",&durmaKosulu);
    point=(function_answer(end)*start-function_answer(start)*end)/(function_answer(end)-function_answer(start));
    if(function_answer(start)*function_answer(end)<0)
    {
        if(durmaKosulu==1)
    {
            while (fabs(function_answer(point)) > tol && i < max_iter && kontrol==1)
        {
            i++;
            point=(function_answer(end)*start-function_answer(start)*end)/(function_answer(end)-function_answer(start));
            printf("start     : %lf\n",start);
            printf("end       : %lf\n",end);
            printf("point     : %lf\n",point);
            printf("f(start)  : %lf\n",function_answer(start));
            printf("f(end)    : %lf\n",function_answer(end));
            printf("f(point)  : %lf\n",function_answer(point));
            printf("iteration : %d\n",i);
            if (function_answer(point) == 0) {
                kontrol=0;
            }
            else if (function_answer(point) * function_answer(start) < 0) {
                end = point;
            }
            else if(function_answer(point) * function_answer(end) < 0)
            {
                start = point;
            }
            else
                kontrol=0;
        }
    }
        i=0;
        if(durmaKosulu==2)
        {
            while ((end - start) / pow(2, i) > tol && i < max_iter && kontrol==1)
        {
            i++;
            point=(function_answer(end)*start-function_answer(start)*end)/function_answer(end)-function_answer(start);
            printf("start     : %lf\n",start);
            printf("end       : %lf\n",end);
            printf("point     : %lf\n",point);
            printf("f(start)  : %lf\n",function_answer(start));
            printf("f(end)    : %lf\n",function_answer(end));
            printf("f(point)  : %lf\n",function_answer(point));
            printf("iteration : %d\n",i);
            if (function_answer(point) == 0) {
                kontrol=0;
            }
            else if (function_answer(point) * function_answer(start) < 0) {
                end = point;
            }
            else if(function_answer(point) * function_answer(end) < 0)
            {
                start = point;
            }
            else
                kontrol=0;
        }
    }
    }
    else
    {
        printf("Bu aralikta kok bulunmamaktadır!");
    }
    printf("Kok: %lf\n", point);
}
void newton_raphson()
{
    double x0,epsilon,temp;
    int maxIter,iter=0,kontrol=1;

    printf("Baslangic tahmini x0: ");
    scanf("%lf", &x0);
    printf("Hata toleransi epsilon: ");
    scanf("%lf", &epsilon);
    printf("Maksimum iterasyon sayisi: ");
    scanf("%d", &maxIter);
    double x = x0;
    while (fabs(function_answer(x)) >= epsilon && iter < maxIter && kontrol==1) {
        if (turev_answer(x) == 0) {
            printf("Sifira bolme hatasi!\n");
            kontrol=0;
        }
        double delta_x = function_answer(x) / turev_answer(x);
        iter++;
        printf("Xn: %lf\n",x + delta_x);
        printf("Xn+1: %lf\n",x);
        printf("f(Xn): %lf\n",function_answer(x));
        printf("f'(Xn): %lf\n",turev_answer(x));
        printf("Iterasyon : %d\n", iter);
        x -= delta_x;
    }

    printf("Sonuc: %lf\n",x);
}
void sayisal_turev()
{
    double x, h, result;
    int secenek;
    printf("x degerini girin: \n");
    scanf("%lf", &x);
    printf("h degerini girin: \n");
    scanf("%lf", &h);
    printf("ileri farklar yontemi icin: 1\n");
    printf("geri farklar yontemi icin: 2\n");
    printf("merkezi farklar yontemi icin: 3\n");
    printf("Secenek: ");
    scanf("%d",&secenek);
    if(secenek==1)
    {
       result = (function_answer(x + h) - function_answer(x)) / h;

    }else if(secenek==2)
    {
        result = (function_answer(x) - function_answer(x - h)) / h;

    }else if(secenek==3)
    {
        result = (function_answer(x + h) - function_answer(x - h)) / (2 * h);

    }else
    {
        printf("Lutfen verilen seceneklerden birini giriniz!");
    }
    printf("result: %lf",result);
}
void simpson()
{
    double x,y,sonuc=0;
   int n,a,b,i,j=1,secenek;
   printf("integral baslangic degerini girin: ");
   scanf("%d",&a);
   printf("integral son degerini girin: ");
   scanf("%d",&b);
   printf("n sayisini girin: ");
   scanf("%d",&n);
   printf("simpson 1/3 icin : 1\n");
   printf("simpson 3/8 icin : 2\n");
   printf("secenek : ");
   scanf("%d",&secenek);
   x = (b-a)/ n;
   if(x<0)
   {
       printf("integral baslangici, sonundan buyuk olmali!");
   }else if(x==0)
   {
       printf("yontem sonucu: 0");
   }else if(secenek==1)
   {
        for(i=0;i<n;i++)
        {
            sonuc += (x * (function_answer(a+x*i) + 4*function_answer((2*a+x*j)/2) + function_answer(a+x*(i+1))) / 6);
            j+=2;
        }
    printf("yontem sonucu: %lf",sonuc);
   }else if(secenek==2)
   {
        j=1;
        y=x/3;
        for(i=0;i<n;i++)
        {
    sonuc += (x * (function_answer(a+x*i) + 3*function_answer(a+y*j) + 3*function_answer(a+y*(j+1)) + function_answer(a+x*(i+1))) / 8);
    j+=3;
        }
    printf("yontem sonucu: %lf",sonuc);
   }

}
void trapez()
{
    double x,sonuc=0;
   int n,a,b,i;
   printf("integral baslangic degerini girin: ");
   scanf("%d",&a);
   printf("integral son degerini girin: ");
   scanf("%d",&b);
   printf("n sayisini girin: ");
   scanf("%d",&n);
   x = (b-a)/ n;

   if(x<0)
   {
       printf("integral baslangici, sonundan buyuk olmali!");
   }else if(x==0)
   {
       printf("yontem sonucu: 0");
   }else
   {
        for(i=0;i<n;i++)
        {
            sonuc += (x * (function_answer(a+x*i) + function_answer(a+x*(i+1))) / 2);
        }
    printf("yontem sonucu: %lf",sonuc);
   }
}
void matris_al()
{
    int i,j;
    printf("Matris boyutunu girin: ");
    scanf("%d",&function1.n);
    for(i=0;i<function1.n;i++)
    {
        for(j=0;j<function1.n;j++)
        {
            printf("[%d][%d]: ",i,j);
            scanf("%lf",&function1.Matris[i][j]);
        }
    }
    printf("Matris: \n");
    for(i=0;i<function1.n;i++)
    {
        for(j=0;j<function1.n;j++)
        {
            printf("  %lf",function1.Matris[i][j]);
        }
        printf("\n");
    }
}
void matris_tersi()
{
    int i,j,x;
    double matris2[function1.n][function1.n],d,k;
    for(i=0;i<function1.n;i++)
    {
        for(j=0;j<function1.n;j++)
        {
            if(i==j)
            {
                matris2[i][j]=1;
            }else
            {
                matris2[i][j]=0;
            }
        }
    }

    for(i=0;i<function1.n;i++)
    {
        d=function1.Matris[i][i];
        for(j=0;j<function1.n;j++)
        {
            if(d!=0)
            {
            function1.Matris[i][j]=function1.Matris[i][j]/d;
            matris2[i][j]=matris2[i][j]/d;
            }
            else
            {
            function1.Matris[i][j]=function1.Matris[i][j];
            matris2[i][j]=matris2[i][j];
            }
        }
        for(x=0;x<function1.n;x++)
        {
            if(x!=i)
            {
                k=function1.Matris[x][i];
                for(j=0;j<function1.n;j++)
                {
                    function1.Matris[x][j]=function1.Matris[x][j]-(function1.Matris[i][j]*k);
                    matris2[x][j] = matris2[x][j] - (matris2[i][j]*k);
                }
            }
        }
    }
    printf("\n matrisin tersi: \n");
    for(i=0;i<function1.n;i++)
    {
        for(j=0;j<function1.n;j++)
        {
            printf(" %lf ",matris2[i][j]);
        }
        printf("\n");
    }

}
void gauss_yoketme()
{
    int i, j, n, k;
    double c, x[diziBoyut], sum = 0.0;
    printf("Matris boyutlarini girin: ");
    scanf("%d", &n);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j <= n; j++)
        {
            printf("[%d][%d]: ", i, j);
            scanf("%lf", &function1.Matris[i][j]);
        }
    }
    printf("Matris: \n");
    for (i = 0; i < n; i++)
    {
        for (j = 0; j <= n; j++)
        {
            printf("  %lf", function1.Matris[i][j]);
        }
        printf("\n");
    }
    for (j = 0; j < n; j++)
    {
        for (i = j + 1; i < n; i++)
        {
            c = function1.Matris[i][j] / function1.Matris[j][j];
            for (k = 0; k <= n; k++)
            {
                function1.Matris[i][k] = function1.Matris[i][k] - c * function1.Matris[j][k];
            }
        }
    }
    x[n - 1] = function1.Matris[n - 1][n] / function1.Matris[n - 1][n - 1];
    for (i = n - 2; i >= 0; i--)
    {
        sum = 0;
        for (j = i + 1; j < n; j++)
        {
            sum = sum + function1.Matris[i][j] * x[j];
        }
        x[i] = (function1.Matris[i][n] - sum) / function1.Matris[i][i];
    }
    printf("\nYontem Sonucu: \n");
    for (i = 0; i < n; i++)
    {
        printf("\nx%d = %lf\t", i + 1, x[i]);
    }
}
void gauss_seidel()
{
    float b[diziBoyut], x[diziBoyut], y[diziBoyut],Matris[diziBoyut][diziBoyut];
    int n = 0, m = 0, i = 0, j = 0;
    printf("NxN'lik matrisin satir sayisini girin: ");
    scanf("%d", &n);
    printf("Iterasyon sayisini girin: ");
    scanf("%d", &m);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("Matrisin [%d][%d]. elemaninin girin: ", i + 1, j + 1);
            scanf("%f", &Matris[i][j]);
        }
    }
    printf("Sonuc matrisinin degerlerini girin:\n");
    for (i = 0; i < n; i++)
    {
        printf("Matrisin [%d]. elemaninin girin: ", i + 1);
        scanf("%f", &b[i]);
    }
    printf("X baslangic degerlerini girin:\n");
    for (i = 0; i < n; i++)
    {
        printf("x%d: ", i + 1);
        scanf("%f", &x[i]);
    }
    while (m > 0)
    {
        for (i = 0; i < n; i++)
        {
            y[i] = (b[i] / Matris[i][i]);
            for (j = 0; j < n; j++)
            {
                if (j == i)
                {
                    continue;
                }
                else
                {
                    y[i] = y[i] - ((Matris[i][j] / Matris[i][i]) * x[j]);
                }
            }
            x[i] = y[i];
            printf("x%d = %f\n", i + 1, y[i]);
        }
        printf("\n");
        m--;
    }
}
void gregory_newton()
{
     double ax[diziBoyut+1], ay [diziBoyut+1], fark[diziBoyut+1][emir+1], nr=1, dr=1,x,p,h,yp;
    int n,i,j,k;
    printf("Girilecek kac adet (x,y) ikilisi bulunmakta : ");
    scanf("%d",&n);


    for (i=0;i<=n;i++)
    {
        printf("%d. x ve y degerlerini girin: ",i+1);
        scanf("%lf %lf",&ax[i],&ay[i]);
    }
    printf("Bulmak istediginiz x degerini giriniz :\n");
    scanf("%lf",&x);
    h=ax[1]-ax[0];
    for (i=0;i<=n-1;i++)
        fark[i][1] = ay[i+1]-ay[i];
    for (j=2;j<=emir;j++)
        for(i=0;i<=n-j;i++)
        fark[i][j] = fark[i+1][j-1] - fark[i][j-1];
    i=0;
    while (!(ax[i]>x))
        i++;
    i--;
    p = (x-ax[i])/h;
    yp = ay[i];
    for (k=1;k<=emir;k++)
    {
        nr =p-k+1;
        dr=k;
        yp +=(nr/dr)*fark[i][k];
    }
    printf("x = %lf icin y = %lf\n",x,yp);
}
int main()
{
    int choice;
    choice = Menu();
    switch(choice)
    {
    case 0:
    break;
    case 1:
    fonk_al();
    bisection();
    break;
    case 2:
    fonk_al();
    regula_falsi();
    break;
    case 3:
    fonk_al();
    turevAl();
    newton_raphson();
    break;
    case 4:
    matris_al();
    matris_tersi();
    break;
    case 5:
    gauss_yoketme();
    break;
    case 6:
    gauss_seidel();
    break;
    case 7:
    fonk_al();
    sayisal_turev();
    break;
    case 8:
    fonk_al();
    simpson();
    break;
    case 9:
    fonk_al();
    trapez();
    break;
    case 10:
    gregory_newton();
    break;
    }
    return 0;
}
