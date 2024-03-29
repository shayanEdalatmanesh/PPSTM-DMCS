
#ifndef IO_cpp
#define IO_cpp

#include <stdio.h>
#include <stdlib.h>
#include <locale.h> 


void str_cp(int n, const char* src, char* dest){
    for(int i=0; i<n; i++){ dest[i]=src[i]; }
    dest[n]='\0';
}


extern "C" {

    int ReadNumsUpTo_C (char *fname, double *numbers, int * dims, int noline) {
        FILE *f;
        char line[5000]; // define a length which is long enough to store a line
        char *waste;
        int waste2;
        long i=0, j=0, k=0, tot=0; 
        int nx=dims[0];
        int ny=dims[1];
        int nz=dims[2];
        printf ("FileRead program: reading %s file\n", fname);
        printf ("XYZ dimensions are %d %d %d\n", dims[0], dims[1], dims[2]);
        f=fopen(fname, "r");
        if (f==NULL)        {
            fprintf(stderr, "Can't open the file %s", fname);
            exit (1); 
        }
        for (i=0; i<noline; i++) {   
            waste=fgets(line,5000, f);
        }
//       printf ("Line: %s", line);
        for  (tot=0, k=0; k<dims[2]; k++){
            for (j=0; j<dims[1]; j++){   
                for (i=0; i<dims[0]; i++){
                    waste2=fscanf(f,"%lf",&numbers[tot]);
//                    printf ("%20.20lf ", numbers[tot]);
                    tot++;
//                    if (tot > 5 ) exit(1);
                }
            }
        }
//       printf ("%lf %lf %lf %lf %lf\n", numbers[tot-1], numbers[tot-2], numbers[tot-3], numbers[tot-4], numbers[tot-5]);
        printf("Reading DONE\n");
        fclose(f);
        return 0;
    }

    int read_AIMS_coefs(char *fname, double* coefs, int* period, int nMOmax, int nMOmin, int nAtoms, int nPerAtoms ){
        setlocale( LC_ALL, "C" ); // https://msdn.microsoft.com/en-us/library/x99tb11d(v=vs.71).aspx
        FILE *f;
        const int nchmax = 1000000;
        char* line = new char[nchmax];
        //char line[nchmax]; // define a length which is long enough to store a line
        char word[64];
        f=fopen(fname, "r");
        int nPerMO = nPerAtoms*nAtoms;

        const int signs[9] = {1, -1,-1,1,   -1,-1,-1,1,-1 }; //# {1, 1, 1, -1, 1, 1, 1, 1, -1, 1};(*Dont change, means - +s, +py +pz -px +dxy +dyz +dz2 -dxz +dx2y2)
                            //# but l=1 has opposite phase than l=0 and l=2 is n-1 - the same phase as l=1 ==>  sign[s]*{1, -1, -1, 1, -1, -1, -1, 1, -1};
        fgets(line,nchmax,f);

        int iline=0;
        while( fgets(line,nchmax,f) ){
            printf( "%s", line );
            iline++;
            if( line[0]=='#' ) continue;
            char ctyp = line[12];  // atomic?
            if( ctyp != 'a' ) continue;
            char l   = line[24];       // l-quantum number
            str_cp(5,line+6 ,word); int ia = strtod(word, NULL)-1; printf( " ia |%s| \n", word); // atom
            str_cp(2,line+21,word); int n  = strtod(word, NULL); printf( " n  |%s| \n", word); // n-quantum number
            str_cp(2,line+27,word); int m  = strtod(word, NULL); printf( " m  |%s| \n", word); // m-quantum number
            int ioff = -1;
            printf( " ia %i n,l,m %i %c %i \n", ia, n, l, m  ); //Debug
            switch(l){
                case 's': if( n==period[ia]   )                     ioff = 0;   break;
                case 'p': if( n==period[ia]   )                     ioff = m+2; break; // -1 -> 1
                case 'd': if((n==period[ia]-1 ) && (nPerAtoms>=9) ) ioff = m+6; break; // -2 -> 4
            }
            if(ioff>=0){
                int sign   = signs[ioff] * ((period[ia] % 2) * 2 - 1); //# phase of radial function in long distance for l=0: if n even - +1, if odd - -1
                printf( "line# %i ia,per %i %i n,l,m %i %c %i io,sign %i %i \n", iline, ia, period[ia], n, l, m, ioff, sign ); //DEBUG1111
                double* coefs_ = coefs + ia*nPerAtoms + ioff;
                char* s = line+34+30*nMOmin;
                for(int imo=0; imo<2*(nMOmax-nMOmin); imo++){
                    str_cp(12,s, word);
                    *coefs_ = strtof(word, NULL) * sign;
                    printf( " (%04i,%04i,%03i) |%s| %g \n", ia, imo, ioff, word, *coefs_  );
                    coefs_ += nPerMO;
                    s+=15;
                }
            }
        }
        fclose(f);
        delete [] line;
        return 0;
    }

}

#endif

