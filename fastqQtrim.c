#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define MAXLETTER 512
#define QVALUEnUM 93
#define HEADtRIM 0
#define TAILtRIM 0
#define LENGTHlIMIT 0
#define QUALITYtRIMMING 25
#define MAXlENGTH 0

typedef struct{
  int headTrim;
  int tailTrim;
  int shortestLength;
  int longestLength;
  int lowestQuality;
  char showQualityHistogram;
}OPTIONS;

typedef struct{
  char *string;
  int buffer;
  int length;
}STRING;

typedef struct{
  char *name;
  char *dna;
  char *quality;
  char *dnaShow;
  char *qualityShow;
  int length;
}FASTQ;

OPTIONS interface(int argc, char *argv[]);
int checkOptionValue(char *key, char *value, int limit);
int readData(OPTIONS options);
int destroyFASTQ(FASTQ *sequence);
FASTQ readOneSequence(FILE *fin);
char *readOneLine(FILE *fin);
int processSequence(OPTIONS options, FASTQ *sequence);
int showResult(FASTQ sequence);
int showHistogram(int *histogram);
int initializeHistogram(int *histogram, int number);
int *cumulateHistogram(int *histogram, FASTQ sequence);
int chomp(char *string);
void help(void);

int main(int argc, char *argv[]){
  OPTIONS options;
  options=interface(argc, argv);
  readData(options);
  return 0;
}

OPTIONS interface(int argc, char *argv[]){
  OPTIONS options;
  int count;
  char *key;
  char *value;
  if((argc % 2) == 0){
    help();
  }

  options.headTrim       = HEADtRIM;
  options.tailTrim       = TAILtRIM;
  options.shortestLength = LENGTHlIMIT;
  options.longestLength  = MAXlENGTH;
  options.lowestQuality  = QUALITYtRIMMING;
  options.showQualityHistogram = 'F';
  for(count=1; count < argc-1; count += 2){
    key=argv[count];
    value=argv[count+1];
    if(strcmp(key, "-h") == 0){
      options.headTrim=checkOptionValue(key, value, 0);
    }else if(strcmp(key, "-t") == 0){
      options.tailTrim=checkOptionValue(key, value, 0);
    }else if(strcmp(key, "-s") == 0){
      options.shortestLength=checkOptionValue(key, value, 0);
    }else if(strcmp(key, "-l") == 0){
      options.longestLength=checkOptionValue(key, value, 0);
    }else if(strcmp(key, "-q") == 0){
      options.lowestQuality=checkOptionValue(key, value, 0);
    }else if(strcmp(key, "-p") == 0){
      if(strcmp(value, "T") != 0 && strcmp(value, "F") != 0){
	fprintf(stderr, "The value for %s should be T/F\n", key);
	help();
      }else{
	options.showQualityHistogram=value[0];
      }
    }else{
      fprintf(stderr, "Unknown option: %s %d\n", key, count);
      help();
    }
  }

  if(options.longestLength > 0 && options.longestLength < options.shortestLength){
    fprintf(stderr, "-s value should not be bigger than -l value\n");
    exit(1);
  }
  return options;
}

int checkOptionValue(char *key, char *value, int limit){
  int number;
  number=atoi(value);
  if(number < limit){
    fprintf(stderr, "Number for %s should be integer >= %d\n"
	   , key, limit);
    help();
  }
  return number;
}

void help(void){
  fprintf(stderr, "\n%s compiled on %s %s\n", __FILE__, __DATE__, __TIME__);

  fprintf(stderr, "  --- usage ---\n");
  fprintf(stderr, " -h: [integer >= 0] Number of bases to trim from the head. \n");
  fprintf(stderr, " -t: [integer >= 0] Number of bases to trim from the tail.\n");
  fprintf(stderr, " -s: [integer >= 0] Length of shortest sequence to show after trimming.\n");
  fprintf(stderr, " -l: [integer >= 0] Sequences longer than this value will be \n");
  fprintf(stderr, "     trimmed to this length. If 0 (zero) trimming will not be performed.\n");
  fprintf(stderr, " -q: [integer >= ] Lowest quality limit. \n");
  fprintf(stderr, "     If this program find a base lower than this value,\n");
  fprintf(stderr, "     the sequence will be trimmed there.\n");
  fprintf(stderr, " -p: [T/F] Just pursue and show quality distribution.\n");
  exit(1);
}

int readData(OPTIONS options){
  FASTQ oneSequence;
  int *histogram;
  if(options.showQualityHistogram=='T'){
    histogram=(int *)calloc(QVALUEnUM, sizeof(int));
  }else{
    histogram=NULL;
  }
  while(!feof(stdin)){
    oneSequence=readOneSequence(stdin);
    if(oneSequence.length==0){
      break;
    }
    if(oneSequence.length!=0 &&
       processSequence(options, &oneSequence)!=0){
      if(histogram!=NULL){
	histogram=cumulateHistogram(histogram, oneSequence);
      }else{
	showResult(oneSequence);
      }
    }else{
      /* fprintf(stderr, "skip\n"); */
    }
    destroyFASTQ(&oneSequence);
  }
  
  if(histogram!=NULL){
    showHistogram(histogram);
    free(histogram);
  }
  return 0;
}

int destroyFASTQ(FASTQ *sequence){
  if(sequence->name != NULL){
    free(sequence->name);
    sequence->name=NULL;
  }
  if(sequence->dna != NULL){
    free(sequence->dna);
    sequence->dna=NULL;
  }
  if(sequence->quality != NULL){
    free(sequence->quality);
    sequence->quality=NULL;
  }
  return 0;
}

FASTQ readOneSequence(FILE *fin){
  FASTQ sequence;
  char *dummy;

  do{
    sequence.name=readOneLine(fin);
  }while(sequence.name[0]!='@' && sequence.name[0]!='\0');
  sequence.dna=readOneLine(fin);
  dummy=readOneLine(fin);
  /*if(strcmp(dummy, "+\n")!=0 && !feof(stdin)){*/
  if(dummy[0]!='+' && !feof(stdin)){
      fprintf(stderr, "This file does not seem to be the fastq format. EXIT\n");
      exit(3);
  }
  free(dummy);
  sequence.quality=readOneLine(fin);
  chomp(sequence.name);
  chomp(sequence.dna);
  chomp(sequence.quality);
  sequence.length=strlen(sequence.dna);
  sequence.dnaShow=NULL;
  sequence.qualityShow=NULL;
  return sequence;
}

char *readOneLine(FILE *fin){
  char tmp[MAXLETTER];
  char *string;
  int tmplength;
  int length=0;
  int buffer=MAXLETTER;

  /* initialize string */
  if((string=(char *)malloc(buffer * sizeof(char)))==NULL){
    fprintf(stderr, "Could not initialize memory on line %d\n"
	    , __LINE__);
    exit(2);
  }
  string[0]='\0';
  
  while((fgets(tmp, MAXLETTER-1, fin))!=NULL){
    tmplength=strlen(tmp);
    /* realloc memory */
    if(tmplength + length > buffer){
      buffer = length + tmplength + MAXLETTER;
      if((string=(char *)realloc(string, buffer * sizeof(char)))==NULL){
	exit(2);
      }
    }

    strcat(string, tmp);
    length+=tmplength;
    if(strstr(tmp, "\n")!=NULL){
      break;
    }
  }
  return string;
}

int processSequence(OPTIONS options, FASTQ *sequence){
  int position;
  if(sequence->length < options.headTrim + options.shortestLength){
    /*
    printf("Too short sequence (%d bp) before trimming. %s\n"
	   , sequence->length, sequence->name);
    */
    return 0;
  }
  /* trim tail */
  position = sequence->length - options.tailTrim;
  sequence->dna[position] = '\0';
  sequence->quality[position] = '\0';
  /* trim tail */
  sequence->dnaShow     = sequence->dna     + options.headTrim * sizeof(char);
  sequence->qualityShow = sequence->quality + options.headTrim * sizeof(char);
  
  for(position=0; sequence->dnaShow[position]!='\0'; position++){
    if(sequence->qualityShow[position] - '!'  < options.lowestQuality){
      sequence->dnaShow[position]='\0';
      sequence->qualityShow[position]='\0';
      break;
    }
  }
  if(options.longestLength > 0 && position > options.longestLength){
    sequence->dnaShow[options.longestLength]='\0';
    sequence->qualityShow[options.longestLength]='\0';
  }
  sequence->length=strlen(sequence->dnaShow);
  if(sequence->length < options.shortestLength){
    /*
    printf("Too short sequence (%d bp) after trimming. %s\n"
	   , sequence->length, sequence->name);
    */
    return 0;
  }else{
    return 1;
  }
}

int showResult(FASTQ sequence){
  printf("%s\n", sequence.name);
  printf("%s\n", sequence.dnaShow);
  printf("+\n");
  printf("%s\n", sequence.qualityShow);
  return 0;
}

int *cumulateHistogram(int *histogram, FASTQ sequence){
  int position;
  int quality;
  for(position=strlen(sequence.qualityShow)-1; position>=0; position--){
    quality=sequence.qualityShow[position] - '!';
    if(quality < 0){
      fprintf(stderr, "something is wrong with quality value\n");
      exit(5);
    }
    histogram[quality]++;
  }
  return histogram;
}

int showHistogram(int *histogram){
  int i;
  for(i=0; i < QVALUEnUM; i++){
    printf("%02d (%c) %d\n", i, i+'!', histogram[i]);
  }
  return 0;
}

int chomp(char *string){
  int position;
  for(position=strlen(string); position >= 0; position--){
    if(iscntrl(string[position])){
      string[position]='\0';
    }else{
      break;
    }
  }
  return 0;
}
