#include "FlagCreator.h"

void WriteFlag(const char* fileName, short** flag, int sizeX, int sizeY)
{
    FILE *datei;
    datei=fopen(fileName,"wb");
    if (datei == NULL)
    {
        printf("Can't open the file.\n");
    }
    if(datei != NULL){
        fwrite(&sizeX,sizeof(int),1,datei);
        fwrite(&sizeY,sizeof(int),1,datei);
        for(int i=0;i<sizeX;i++){
                fwrite(flag[i],sizeof(short),sizeY,datei);
            }
    }
    fclose (datei);
}

short** ReadFlag(const char* fileName, int* sizeX, int* sizeY)
{
    FILE* datei;
    datei=fopen(filename,"rb");
    if(datei == NULL){
        printf("Can't open the file.\n");
    }
    if (datei != NULL){
        fread(sizeX,sizeof(int),1,datei);
        fread(sizeY,sizeof(int),1,datei);
        int** field=mallocFlag(*sizeX,*sizeY);
        for(int i=0;i<*sizeX;i++){
            fread(field[i],sizeof(short),*sizeY,datei);
        }
    }
    fclose(datei);
    return field;
}

short** mallocFlag(int sizeX, int sizeY)
{
    short** field=malloc((sizeX)*sizeof(short*));
    for (int i=0; i<sizeX;i++){
        field[i]=malloc((sizeY)*sizeof(short));
    }
    return field;
}

void CreateFlag()
{
    int sizeX;
    int sizeY;
    printf("\nFirst you need to create a field. Enter size in x-direction: ");
    do {scanf("%d",&sizeX);} while ( getchar() != '\n' );
    printf("\nNow enter size in y-direction: ");
    do {scanf("%d",&sizeY);} while ( getchar() != '\n' );
    short** field=mallocFlag(sizeX,sizeY);
    char run;
    char answer;
    int isize;
    int jsize;
    int ipos;
    int jpos;
    short correct;
    char* filepath;
    double xlength;
    double ylength;
    printf("\nYou created a field. Would you like to add an obstacle? Enter Y or N: ");
    do {scanf("%c",&run);} while ( getchar() != '\n' );
    if(run == Y){
        printf("\nThe field has size %d*%d.Please remember the rules for obstacles!",sizeX,sizeY);
    }
    printf("\nWe need more Information for Visualization. Please enter xlength: ");
    do {scanf("%lf",&xlength);} while ( getchar() != '\n' );
    printf("\nAnd now ylength: ");
    do {scanf("%lf",&ylength);} while ( getchar() != '\n' );
    double delx=xlength/(double)sizeX;
    double dely=ylength/(double)sizeY;
    while(run == Y){
        printf("\nWhich size should the obstacle have in x-direction? Enter a number: ");
        do {scanf("%d",&isize);} while ( getchar() != '\n' );
        while(isize>sizeX){
            printf("\nPick a smaller size please. Enter a new number: ");
            do {scanf("%d",&isize);} while ( getchar() != '\n' );
        }
        printf("\nWhich size should the obstacle have in y-direction? Enter a number: ");
        do {scanf("%d",&jsize);} while ( getchar() != '\n' );
        while(jsize>sizeY){
            printf("\nPick a smaller size please. Enter a new number: ");
            do {scanf("%d",&jsize);} while ( getchar() != '\n' );
        }
        printf("\nSet the position for the lower left block of the obstacle: in x-direction: ");
        do {scanf("%d",&ipos);} while ( getchar() != '\n' );
        printf("\nAnd now in y-direction:\ ");
        do {scanf("%d",&jpos);} while ( getchar() != '\n' );
        if(ipos<1 | jpos<1 | (ipos+size)>sizeX | (jpos+size) > sizeY){
            printf("\nWrong positions or size too big. No obstacle added.");
        }else{
            for(int i=ipos;i<=ipos+size;i++){
                for(int j=jpos;j<=jpos+size;j++){
                    field[i-1][j-1]=1;
                }
            }
        }
        printf("\nDo you want to check your field in Paraview? Enter Y or N: ");
        do {scanf("%c",&answer);} while ( getchar() != '\n' );
        if(answer == Y){
            printf("\nEnter a filepath to save the file:");
            do {scanf("%s",filepath);} while ( getchar() != '\n' );
            writeVTKfileFor2DintegerField(filepath,"shows obstacles",field,sizeX,sizeY,delx,dely);
        }
        printf("\nDo you want to add another obstacle? Enter Y or N: ");
        do {scanf("%c",&run);} while ( getchar() != '\n' );
    }
    printf("\nYou have created a field with obstacles. Do you wanna check if it's correct? Enter Y or N: ");
    do {scanf("%c",&answer);} while ( getchar() != '\n' );
    if(answer == Y){
            correct=CorrectnessCheck(field,sizeX,sizeY);
    }
    if(correct == 0){
        printf("\nYour field is not correct. Check the error(s) above and edit it.");
        EditFlag(field,sizeX,sizeY);
    }
    if(correct == 1){
        printf("\nYour field is correct!");
        }
    printf("\nDo you want to save your field? Enter Y or N: ");
    do {scanf("%c",&answer);} while ( getchar() != '\n' );
    if(answer == Y){
        printf("\nEnter a filepath to save the file:");
        do {scanf("%s",filepath);} while ( getchar() != '\n' );
        WriteFlag(filepath,field,sizeX,sizeY);
    }
    printf("\nDo you want to save your field as .vtk File?Enter Y or N: ");
    do {scanf("%c",&answer);} while ( getchar() != '\n' );
    if(answer == Y){
        printf("\nEnter a filepath to save the file:");
        do {scanf("%s",filepath);} while ( getchar() != '\n' );
        writeVTKfileFor2DintegerField(filepath,"shows obstacles",field,sizeX,sizeY,delx,dely);
    }
    printf("\nYou're finished!");
}

short CorrectnessCheck(short** flag, int sizeX, int sizeY)
{

}

void EditFlag(short** flag, int sizeX, int sizeY)
{

}
