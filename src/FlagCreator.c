#include "FlagCreator.h"

void WriteFlag(const char* fileName, int** flag, int sizeX, int sizeY)
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
                fwrite(flag[i],sizeof(int),sizeY,datei);
            }
    }
    fclose (datei);
}

int** ReadFlag(const char* fileName, int* sizeX, int* sizeY)
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
            fread(field[i],sizeof(int),*sizeY,datei);
        }
    }
    fclose(datei);
    return field;
}

int** mallocFlag(int sizeX, int sizeY)
{
    int** field=malloc((sizeX)*sizeof(int*));
    for (int i=0; i<sizeX;i++){
        field[i]=malloc((sizeY)*sizeof(int));
    }
    return field;
}

int** CreateFlag(int sizeX, int sizeY)
{
    int** field=mallocFlag(sizeX,sizeY);
    char run;
    char answer;
    int size;
    int ipos;
    int jpos;
    int correct=1;
    char* filepath;
    printf("\nYou created a field. Would you like to add an obstacle? Enter Y or N: ");
    do {scanf("%c",&run);} while ( getchar() != '\n' );
    if(run == Y){
        printf("\nThe field has size %d*%d.Please remember the rules for obstacles!",sizeX,sizeY);
    }
    while(run == Y){
        printf("\nWhich size should the square have?Enter a number: ");
        do {scanf("%d",&size);} while ( getchar() != '\n' );
        if(size>sizeX | size>sizeY){
            printf("\nPick a smaller size please. Enter a new number: ");
            do {scanf("%d",&size);} while ( getchar() != '\n' );
        }
        printf("\nSet the position for the lower left block of the obstacle: in x-direction: ");
        do {scanf("%d",&ipos);} while ( getchar() != '\n' );
        printf("\nAnd now in y-direction:\ ");
        do {scanf("%d",&jpos);} while ( getchar() != '\n' );

        //fehlermeldung FEHLT
        //obstacle einfügen FEHLT

        printf("\nYou added an obstacle! Do you want to add another?Enter Y or N: ");
        do {scanf("%c",&run);} while ( getchar() != '\n' );
    }
    printf("\nYou have created a field with obstacles. Do you wanna check if it's correct?Enter Y or N: ");
    do {scanf("%c",&answer);} while ( getchar() != '\n' );
    if(answer == Y){

            //korrektheischeck FEHLT

    }
    if(correct == 0){
        printf("\nYour field is not correct. Check the error(s) above and edit it.");

        // FIELD edit FEHLT

    }
    if(correct == 1){
        printf("\nYour field is correct!");
        }
    printf("\nDo you want to save your field?Enter Y or N: ");
    do {scanf("%c",&answer);} while ( getchar() != '\n' );
    if(answer == Y){
        printf("\nEnter a filepath to save the file:");
        do {scanf("%s",filepath);} while ( getchar() != '\n' );
        WriteFlag(filepath,field,sizeX,sizeY);
    }
}
