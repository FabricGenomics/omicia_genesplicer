/* 
 *  Copyright (c) 2003, The Institute for Genomic Research (TIGR), Rockville,
 *  Maryland, U.S.A.  All rights reserved.

 *   genesplicer.cpp was designed by Mihaela PERTEA to find splice 
 *   sites in a fasta file
 
 * Modified by Omicia inc. to read sequences from STDIN as well
 

*/

#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <ctype.h>
#include <cstring>
#include <iostream>
#include <vector>
#ifndef PATH_MAX
#define PATH_MAX 255
#endif


// the following values are directly imported from human/config_file
double ACCEPTOR_HIGH_LIMIT = 14.162281;
double DONOR_HIGH_LIMIT = 13.220582;;
double ACCEPTOR_THRESHOLD = 1.016684;
double DONOR_THRESHOLD = 0.463981;
int DISTA = 100;  // represents the distance within I choose the best acceptor score
int DISTD = 150;  // represents the distance within I choose the best donor score
unsigned char donor_tree = '1';
unsigned char acceptor_tree = '1';

char TRAIN_DIR[PATH_MAX+1]="";

unsigned char use_cod_noncod;

struct Site {
    int type;
    long pos;
    double score;
    Site *prev;
    Site *forw;
    int dir;
};

struct output_struct {
    Site *List;
    Site *Cap;
};

//enum SpliceType { ACC, DON, FALSE};
//enum Direction {F,R};

char FILE_NAME[100];
int UseFile=0;
bool use_stdin = false;


FILE * File_Open  (const char *, const char *);
void LoadData(FILE *,char *);
void  Process_Options  (int, char * []);
output_struct AnalyzeData(char *);



int main(int argc, char * argv []) {
    FILE  * fp, *ofp;
    char  Data[1000001];
    std::string sData;
    long int  Data_Len;
    Site *List,*nod, *site, *Cap;
    void printsites(Site *, Site *, FILE *);
    int B[200],k;
    output_struct output;
    size_t pos = 0;

    std::vector<std::string> vData;

    long offset,i,j;
    char filename[PATH_MAX+1];
    float readdouble;

    if  (argc < 3 ) {
        fprintf (stderr,
                 "USAGE:  %s <fasta-file | STDIN> <specific-genome-training-directory> [options] \n",
                 argv [0]);
        if(argc==2) {
            if(argv[1][0]=='-' && argv[1][1]=='h') {
                printf("Options:\n");
                printf("-f file_name     Write the results in file_name\n");
                printf("-a t             Choose t as a threshold for the acceptor sites\n");
                printf("-d t             Choose t as a threshold for the donor sites\n");
                printf("-e n             The maximum acceptor score within n bp is chosen\n");
                printf("-i n             The maximum donor score within n bp is chosen\n");
                printf("-h               Display the options of the program\n");
            }
        }
        exit (0);
    }
    strcpy(TRAIN_DIR,argv[2]);

     // read config file
//    sprintf(filename,"%s/config_file",TRAIN_DIR);
//    fp = fopen (filename, "r");
//    if  (fp == NULL) {
//        fprintf (stderr, "ERROR:  Could not open config file  %s \n", filename);
//        exit (EXIT_FAILURE);
//    }


//    fscanf(fp,"%f",&readdouble);ACCEPTOR_HIGH_LIMIT=(double)readdouble;
//    fscanf(fp,"%f",&readdouble);DONOR_HIGH_LIMIT=(double)readdouble;
//    fscanf(fp,"%f",&readdouble);ACCEPTOR_THRESHOLD=(double)readdouble;
//    fscanf(fp,"%f",&readdouble);DONOR_THRESHOLD=(double)readdouble;
//    fscanf(fp,"%d",&k); donor_tree=(unsigned char)k;
//    fscanf(fp,"%d",&k); acceptor_tree=(unsigned char)k;
//    fscanf(fp,"%d",&DISTA);
//    fscanf(fp,"%d",&DISTD);
    //  fscanf(fp,"%u",&use_cod_noncod);
//    fclose(fp);

    Process_Options (argc, argv);   // Set global variables to reflect status of
                                  // command-line options.

    if (std::string(argv[1]) == "STDIN" or std::string(argv[1]) == "stdin" )
        use_stdin = true;
    else
        fp = File_Open (argv[1], "r");

/*    cout << "Usefile " << usefile << endl;*/
    if(UseFile) {
        ofp= File_Open (FILE_NAME, "w");
    }

    offset=1;

    
    if (!use_stdin) {
        while (!feof(fp)) {
            LoadData(fp, Data);
            output = AnalyzeData(Data);
            // print splice sites found so far
            if (UseFile)
                printsites(output.List, output.Cap, ofp);
            else
                printsites(output.List, output.Cap, stdout);
            offset += Data_Len-161;
            fprintf(stderr,"Done %ldbp..............\n", offset+160);
        }
    }
    else {
//        std::cin.ignore(26, '\n'); //ignore 26 characters or to a newline, whichever comes first
        std::cin >> sData;
        std::string delimiter = "!";
        while ((pos = sData.find(delimiter)) != std::string::npos) {
            vData.push_back(sData.substr(0, pos));
            sData.erase(0, pos + delimiter.length());
        }
        vData.push_back(sData);

        for (int i=0; i < vData.size(); i++) {
            if (vData[i][0] == '>') {
                fprintf(stdout, "%s\n", vData[i].c_str());
            }
            else {
                std::strcpy(Data, vData[i].c_str());
                output = AnalyzeData(Data);
                printsites(output.List, output.Cap, stdout);
            }
        }
    }

    if (!use_stdin)
      fclose(fp);
    if (UseFile)
      fclose(ofp);
}


output_struct AnalyzeData(char *Data) {

    Site *List,*nod, *site, *Cap;
    char  Filter(char );

    int  Is_Acceptor  (const int *, double *);
    int  Is_Donor  (const int *, double *);

    double score;
    long int offset,i,j;
    long int Data_Len;
    int B[200],k;
    
    List=NULL;
    nod=List;
    Cap=List;
    output_struct output;
    offset=1;
    
    Data_Len = strlen (Data);
    
    //for debugging purposes
    //    printf("%ldbp\n", Data_Len);
    //printf(Data);
    //printf("\n");
    
    
    for  (i = 0;  i < Data_Len;  i ++) {
        // Converts all characters to  acgt
        Data [i] = Filter (Data [i]);
    }
    
    for  (i = 80;  i <= Data_Len-82;  i ++) {
        // look fot gt and ag
        
        // forward direction
        if (Data[i]=='a' && Data[i+1]=='g') { // Deal with acceptors
            k=0;
            for (j=i-80;j<i+82;j++) {
                switch (Data[j]){
                    case 'a': B[k]=0;break;
                    case 'c': B[k]=1;break;
                    case 'g': B[k]=2;break;
                    case 't': B[k]=3;break;
                    default: B[k]=1;
                }
                k++;
            }
            
            if (Is_Acceptor(B,&score)) {
                //	    printf("ag: %ld\n",i+offset);
                
                // add acceptor to list;
                site=(Site *) malloc(sizeof(Site));
                if (site == NULL) {
                    fprintf(stderr,"Memory allocation for Site position failure.\n");
                    abort();
                }

                site->type=1;
                site->score=score;
                site->pos=i+offset;
                site->dir=1;
                site->forw=NULL;
                site->prev=nod;
                
                if (List == NULL)
                    List=site;
                Cap = site;
                if (nod != NULL )
                    nod->forw = site;
                nod = site;
            }
        }

        if (Data[i]=='g' && Data[i+1]=='t') { // Deal with donors
            k=0;
            for(j=i-80;j<i+82;j++){
                switch (Data[j]){
                    case 'a': B[k]=0;break;
                    case 'c': B[k]=1;break;
                    case 'g': B[k]=2;break;
                    case 't': B[k]=3;break;
                    default: B[k]=1;
                }
                k++;
            }
            
            if(Is_Donor(B,&score)) {
                //printf("gt: %ld\n",i+offset+1);fflush(stdout);
                
                // add donor to list;
                site=(Site *) malloc(sizeof(Site));
                if (site == NULL) {
                    fprintf(stderr,"Memory allocation for Site position failure.\n");
                    abort();
                }
                
                site->type=2;
                site->score=score;
                site->pos=i+offset;
                site->dir=1;
                site->forw=NULL;
                site->prev=nod;
                
                if (List == NULL)
                    List = site;
                Cap=site;
                if (nod != NULL)
                    nod->forw=site;
                nod=site;
                
            }
        }


        /*
         * Omicia's input sequences are always forward stranded.
         * Analyzing the revard strand is not necessary
         * changed by bstade
         */

        // reversed direction
//        if(Data[i]=='c' && Data[i+1]=='t') { // Deal with acceptors
//            k=0;
//            for(j=i+81;j>=i-80;j--){
//                switch (Data[j]){
//                    case 'a': B[k]=3;break;
//                    case 'c': B[k]=2;break;
//                    case 'g': B[k]=1;break;
//                    case 't': B[k]=0;break;
//                    default: B[k]=1;
//                }
//                k++;
//            }
//
//            if(Is_Acceptor(B,&score)) {
//                //	    printf("ag: %ld\n",i+offset+1);
//
//                // add acceptor to list;
//                site=(Site *) malloc(sizeof(Site));
//                if (site == NULL) {
//                    fprintf(stderr,"Memory allocation for Site position failure.\n");
//                    abort();
//                }
//
//                site->type=1;
//                site->score=score;
//                site->pos=i+offset;
//                site->dir=-1;
//                site->forw=NULL;
//                site->prev=nod;
//
//                if (List == NULL) List=site;
//                Cap=site;
//                if (nod != NULL) nod->forw=site;
//
//                nod=site;
//
//            }
//        }
//
//        if (Data[i]=='a' && Data[i+1]=='c') { // Deal with donors
//
//            k=0;
//            for(j=i+81;j>=i-80;j--){
//                switch (Data[j]){
//                    case 'a': B[k]=3;break;
//                    case 'c': B[k]=2;break;
//                    case 'g': B[k]=1;break;
//                    case 't': B[k]=0;break;
//                    default: B[k]=1;
//                }
//                k++;
//            }
//
//            if(Is_Donor(B,&score)) {
//
//                //	    printf("gt: %ld\n",i+offset+1);
//
//                // add donor to list;
//                site=(Site *) malloc(sizeof(Site));
//                if (site == NULL) {
//                    fprintf(stderr,"Memory allocation for Site position failure.\n");
//                    abort();
//                }
//
//                site->type=2;
//                site->score=score;
//                site->pos=i+offset;
//                site->dir=-1;
//                site->forw=NULL;
//                site->prev=nod;
//
//                if(List == NULL) List=site;
//                Cap=site;
//                if(nod !=NULL) nod->forw=site;
//
//                nod=site;
//            }
//
//        }
        
    }
    output.List = List;
    output.Cap= Cap;
    return output;
}


void LoadData(FILE *fp, char *Data) {
    char ch;
    long length,lcopy,back;
    char line[1001];
    bool newlineFlag;
    newlineFlag = false;

    // find the start
    length=161;
    while(length) {
        fseek(fp,-1,SEEK_CUR);
        ch=fgetc(fp);
        fseek(fp,-1,SEEK_CUR);
        if(ch != '\n') length--;
    }

    // copy the data
    length=0;
    Data[0]='\0';


    while(length<1000000) {
        if(!fgets(line,1000,fp)) break;

        lcopy=strlen(line);

        newlineFlag = false;
        if(line[lcopy-1]=='\n') {
            line[lcopy-1]='\0';
            lcopy = lcopy - 1;
            newlineFlag = true;
        }
        if(length+lcopy>1000000) {
            back=lcopy;
            lcopy=1000000-length;
            back-=lcopy;
            while(back) {
                fseek(fp,-1,SEEK_CUR);
                ch=fgetc(fp);
                fseek(fp,-1,SEEK_CUR);
                if(ch != '\n') back--;
            }
        }
        //if there's a newline, need to subtract one from length and from the lcopy when strncatting; bug fixed by Roger Craig (University of Delaware)
        if(newlineFlag){
            length+=lcopy;
            length-=1;
            strncat(Data,line,lcopy-1);
        }
        else {
            length+=lcopy;
            strncat(Data,line,lcopy);
        }
    }
    Data[length]='\0';
}


void printsites(Site *List, Site *Cap,FILE *ofp)
{
  char confidence[15];
  Site *nod, *scan;

  nod=List;

  while(nod!=NULL) {
    scan=nod->prev;

    if(nod->type==1) {
      while(scan && nod->pos-scan->pos<=DISTA) {
	if(scan->type==1 && scan->score<nod->score) {
	  scan->type=0;
	}
	scan=scan->prev;
      }
    }
      
    if(nod->type==2) {
      while(scan && nod->pos-scan->pos<=DISTD) {
	if(scan->type==2 && scan->score<nod->score) {
	  scan->type=0;
	}
	scan=scan->prev;
      }
    }

    nod=nod->forw;
  }

  nod=Cap;

  while(nod!=NULL) {
    scan=nod->forw;
      
    if(nod->type==1) {
      while(scan && scan->pos-nod->pos<=DISTA) {
	if(scan->type==1 && scan->score<nod->score) {
	  scan->type=0;
	}
	scan=scan->forw;
      }
    }  

    if(nod->type==2) {
      while(scan && scan->pos-nod->pos<=DISTD) {
	if(scan->type==2 && scan->score<nod->score) {
	  scan->type=0;
	}
	scan=scan->forw;
      }
    }
      
    nod=nod->prev;
  }
  

  nod=List;

  while(nod!=NULL) {
    if(nod->type == 1) {
      if(nod->score>=ACCEPTOR_HIGH_LIMIT) strcpy(confidence,"High");
      else strcpy(confidence,"Medium");
      fprintf(ofp,"%ld %d %.6f %s acceptor\n", nod->pos, nod->dir, nod->score, confidence);
    }
    if(nod->type == 2) {
      if(nod->score>=DONOR_HIGH_LIMIT) strcpy(confidence,"High");
      else strcpy(confidence,"Medium");
      fprintf(ofp,"%ld %d %.6f %s donor\n",nod->pos, nod->dir, nod->score, confidence);
    }
    List=nod;
    nod=nod->forw;
    free(List);
  }
}

char  Filter  (char Ch)

//  Return a single  a, c, g or t  for  Ch .  Choice is to minimize likelihood
//  of a stop codon on the primary strand.

{
  switch  (tolower (Ch))
    {
    case  'a' :
    case  'c' :
    case  'g' :
    case  't' :
      return  tolower(Ch);
    case  'r' :     // a or g
      return  'g';
    case  'y' :     // c or t
      return  'c';
    case  's' :     // c or g
      return  'c';
    case  'w' :     // a or t
      return  't';
    case  'm' :     // a or c
      return  'c';
    case  'k' :     // g or t
      return  't';
    case  'b' :     // c, g or t
      return  'c';
    case  'd' :     // a, g or t
      return  'g';
    case  'h' :     // a, c or t
      return  'c';
    case  'v' :     // a, c or g
      return  'c';
    default :       // anything
      return  'c';
    }
}

FILE *  File_Open  (const char * Filename, const char * Mode)

/* Open  Filename  in  Mode  and return a pointer to its control
*  block.  If fail, print a message and exit. */
{
  FILE  *  fp;
  char line[1000];
  int length;
  
  fp = fopen (Filename, Mode);
  if  (fp == NULL)
    {
      fprintf (stderr, "ERROR:  Could not open file  %s \n", Filename);
      exit (1);
    }

  if(Mode[0] == 'r') {
    fgets(line,1000,fp);
  
    if(line[0] != '>') {
      fprintf (stderr,"ERROR: File %s is not a fasta file\n",Filename);
      exit (1);
    }

    length=0;

    while(length<161) {
      if(fgetc(fp) != '\n') length++;
    }
    
  }
  
  return  fp;
}


void  Process_Options  (int argc, char * argv [])

//  Process command-line options and set corresponding global switches
//  and parameters.
//
//    -f file_name     Write results in file_name
//    -d dist          The maximum score within dist bp is chosen          
//    -h               Display the help of the program
{
    int i;

    for  (i = 3;  i < argc;  i ++) {
        switch  (argv [i] [0]) {
        case  '-' :
            switch  (argv [i] [1]) {
            case  'f':
                if (i+1 == argc) {
                    fprintf (stderr, "After -f you must specify output file.\n");
                    exit(1);
                }
                else {
                    strcpy(FILE_NAME,argv[++i]);
                    UseFile=1;
                }
                break;
            case 'a':
                if (i+1 == argc) {
                    fprintf (stderr, "You must specify a threshold after -a.\n");
                    exit(1);
                }
                else {
                    ACCEPTOR_THRESHOLD=atof(argv[++i]);
                }
                break;
            case 'd':
                if (i+1 == argc) {
                    fprintf (stderr, "You must specify a threshold after -d.\n");
                    exit(1);
                }
                else {
                    DONOR_THRESHOLD=atof(argv[++i]);
                }
                break;
            case 'e':
                if(i+1 == argc){
                    fprintf (stderr, "You must specify a distance after -e.\n");
                    exit(1);
                }
                else {
                    DISTA=atoi(argv[++i]);
                }
                break;
            case 'i':
                if (i+1 == argc) {
                    fprintf (stderr, "You must specify a distance after -i.\n");
                    exit(1);
                }
                else {
                 DISTD=atoi(argv[++i]);
                }
                break;
            case 'h':
                printf("Options:\n");
                printf("-f file_name     Write the results in file_name\n");
                printf("-a t             Choose t as a threshold for the acceptor sites\n");
                printf("-d t             Choose t as a threshold for the donor sites\n");
                printf("-e n             The maximum acceptor score within n bp is chosen\n");
                printf("-i n             The maximum donor score within n bp is chosen\n");
                printf("-h               Display the options of the program\n");
                break;
            default:
               fprintf (stderr, "Unrecognized option %s\n", argv [i]);
            }
	    break;
        default:
           fprintf (stderr, "Unrecognized option %s\n", argv [i]);
        }
    }
}

