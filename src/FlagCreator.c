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
        if(ipos<1 | jpos<1 | (ipos+isize)>sizeX | (jpos+jsize) > sizeY){
            printf("\nWrong positions or size too big. No obstacle added.");
        }else{
            for(int i=ipos;i<ipos+isize;i++){
                for(int j=jpos;j<jpos+jsize;j++){
                    field[i-1][j-1]=C_B;
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
    while(correct == 1){
        printf("\nYour field is not correct. Check the error(s) above and edit it.");
        EditFlag(field,sizeX,sizeY,delx,dely);
        correct=CorrectnessCheck(field,sizeX,sizeY);
    }
    if(correct == 0 && answer == Y){
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
    for(int i=2;i<sizeX;i++){
        for(int j=2;j<sizeY;j++){
            if(flag[i][j]==C_B){
                if(flag[i+1][j]==C_F && flag[i-1][j]==C_F){
                    printf("\nError: (%d,%d) has fluid to the east and west.",i,j);
                    return 1;
                }
                if(flag[i][j+1]==C_F && flag[i][j-1]==C_F){
                    printf("\nError: (%d,%d) has fluid to the north and south.",i,j);
                    return 1;
                }
            }
        }
    }
    return 0;
}

void EditFlag(short** flag, int sizeX, int sizeY, double dx, double dy)
{
    char run=Y;
    char answer;
    int isize;
    int jsize;
    int ipos;
    int jpos;
    char* filepath;
    printf("\nYou're editing an existing field. You can remove and add obstacles.");
    while(run == Y){
        printf("\nDo you want to save the field as .vtk for visualization? Enter Y or N");
        do {scanf("%c",&answer);} while ( getchar() != '\n' );
        if(answer == Y){
            printf("\nEnter a filepath to save the file:");
            do {scanf("%s",filepath);} while ( getchar() != '\n' );
            writeVTKfileFor2DintegerField(filepath,"shows obstacles",flag,sizeX,sizeY,dx,dy);
        }
        printf("\nDo you want to add or remove an obstacle? Enter A or R: ");
        do {scanf("%c",&answer);} while ( getchar() != '\n' );
        if(answer == A){
            printf("\nWhich size should the new obstacle have in x-direction? Enter a number: ");
            do {scanf("%d",&isize);} while ( getchar() != '\n' );
            while(isize>sizeX){
                printf("\nPick a smaller size please. Enter a new number: ");
                do {scanf("%d",&isize);} while ( getchar() != '\n' );
            }
            printf("\nWhich size should the new obstacle have in y-direction? Enter a number: ");
            do {scanf("%d",&jsize);} while ( getchar() != '\n' );
            while(jsize>sizeY){
                printf("\nPick a smaller size please. Enter a new number: ");
                do {scanf("%d",&jsize);} while ( getchar() != '\n' );
            }
            printf("\nSet the position for the lower left block of the new obstacle: in x-direction: ");
            do {scanf("%d",&ipos);} while ( getchar() != '\n' );
            printf("\nAnd now in y-direction:\ ");
            do {scanf("%d",&jpos);} while ( getchar() != '\n' );
            if(ipos<1 | jpos<1 | (ipos+isize)>sizeX | (jpos+jsize) > sizeY){
                printf("\nWrong positions or size too big. No obstacle added.");
            }else{
                for(int i=ipos;i<ipos+size;i++){
                    for(int j=jpos;j<jpos+size;j++){
                        flag[i-1][j-1]=C_B;
                    }
                }
            }
        }
        if(answer == R){
            printf("\nOf which size do you want to remove in x-direction? Enter a number: ");
            do {scanf("%d",&isize);} while ( getchar() != '\n' );
            while(isize>sizeX){
                printf("\nPick a smaller size please. Enter a new number: ");
                do {scanf("%d",&isize);} while ( getchar() != '\n' );
            }
            printf("\nOf which size do you want to remove in y-direction? Enter a number: ");
            do {scanf("%d",&jsize);} while ( getchar() != '\n' );
            while(jsize>sizeY){
                printf("\nPick a smaller size please. Enter a new number: ");
                do {scanf("%d",&jsize);} while ( getchar() != '\n' );
            }
            printf("\nSet the position for the lower left block of the removal: in x-direction: ");
            do {scanf("%d",&ipos);} while ( getchar() != '\n' );
            printf("\nAnd now in y-direction:\ ");
            do {scanf("%d",&jpos);} while ( getchar() != '\n' );
            if(ipos<1 | jpos<1 | (ipos+isize)>sizeX | (jpos+jsize) > sizeY){
                printf("\nWrong positions or size too big. No obstacle removed.");
            }else{
                for(int i=ipos;i<ipos+isize;i++){
                    for(int j=jpos;j<jpos+jsize;j++){
                        flag[i-1][j-1]=C_F;
                    }
                }
            }
        }
        printf("Do you want to continue editing? Enter Y or N: ");
        do {scanf("%c",&run);} while ( getchar() != '\n' );
    }
    printf("\nDo you want to save your field? Enter Y or N: ");
    do {scanf("%c",&answer);} while ( getchar() != '\n' );
    if(answer == Y){
        printf("\nEnter a filepath to save the file:");
        do {scanf("%s",filepath);} while ( getchar() != '\n' );
        WriteFlag(filepath,flag,sizeX,sizeY);
    }
}
