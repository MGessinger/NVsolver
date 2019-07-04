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
    datei=fopen(fileName,"rb");
    short **field;
    if(datei == NULL){
        printf("Can't open the file.\n");
    }
    if (datei != NULL){
        fread(sizeX,sizeof(int),1,datei);
        fread(sizeY,sizeof(int),1,datei);
        field=mallocFlag(*sizeX,*sizeY);
        for(int i=0;i<*sizeX;i++){
            fread(field[i],sizeof(short),*sizeY,datei);
        }
    }
    fclose(datei);
    return field;
}

short** mallocFlag(int sizeX, int sizeY)
{
    short** field;
    field=malloc((sizeX)*sizeof(short*));
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
    char run='Y';
    char answer;
    int isize;
    int jsize;
    int ipos;
    int jpos;
    short correct=0;
    char* filepath;
    double xlength;
    double ylength;
    char obstacle;
    printf("\nYou created a field. Would you like to add an obstacle? Enter Y or N: ");
    do {scanf("%c",&run);} while ( getchar() != '\n' );
    if((run == 'Y') || (run == 'y')){
        printf("\nThe field has size %d*%d.Please remember the rules for obstacles!",sizeX,sizeY);
    }
    printf("\nWe need more Information for Visualization. Please enter xlength: ");
    do {scanf("%lf",&xlength);} while ( getchar() != '\n' );
    printf("\nAnd now ylength: ");
    do {scanf("%lf",&ylength);} while ( getchar() != '\n' );
    double delx=xlength/(double)sizeX;
    double dely=ylength/(double)sizeY;
    while((run == 'Y') || (run == 'y')){
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
        printf("\nWhat type of obstacle do you want to add? Press h for help: ");
        do {scanf("%c",&obstacle);} while ( getchar() != '\n' );
        if(obstacle=='h'|| obstacle == 'H'){
            printf("\nYou can choose between the following types of obstacles:\nR for rectangle\nC for Creative Mode");
            printf("\nPick a type:");
            do {scanf("%c",&obstacle);} while ( getchar() != '\n' );
        }

        switch(obstacle)
        {
        case 'R': case 'r':
            printf("\nSet the position for the lower left block of the rectangle: in x-direction: ");
            do {scanf("%d",&ipos);} while ( getchar() != '\n' );
            printf("\nAnd now in y-direction: ");
            do {scanf("%d",&jpos);} while ( getchar() != '\n' );
            if((ipos<1) | (jpos<1) | ((ipos+isize)>sizeX) | ((jpos+jsize) > sizeY)){
                printf("\nWrong positions or size too big. No rectangle added.");
            }else{
                AddRectangle(field,isize,jsize,ipos,jpos,C_B);
            }
            break;

        default:
            printf("\nNo valid input.");
            break;
        }
        printf("\nDo you want to check your field in Paraview? Enter Y or N: ");
        do {scanf("%c",&answer);} while ( getchar() != '\n' );
        if(answer == 'Y' || answer == 'y'){
            printf("\nEnter a filepath to save the file:");
            do {scanf("%s",filepath);} while ( getchar() != '\n' );
            writeVTKfileFor2DintegerField(filepath,"shows obstacles",field,sizeX,sizeY,delx,dely);
        }
        printf("\nDo you want to add another obstacle? Enter Y or N: ");
        do {scanf("%c",&run);} while ( getchar() != '\n' );
    }
    printf("\nYou have created a field with obstacles. Do you wanna check if it's correct? Enter Y or N: ");
    do {scanf("%c",&answer);} while ( getchar() != '\n' );
    if((answer == 'Y') || (answer == 'y')){
            correct=CorrectnessCheck(field,sizeX,sizeY);
    }
    while(correct > 0){
        printf("\nYour field is not correct. Check the error(s) above and edit it.");
        EditFlag(field,sizeX,sizeY,delx,dely);
        correct=CorrectnessCheck(field,sizeX,sizeY);
    }
    if((correct == 0 && answer == 'Y')||(correct == 0 && answer == 'y')){
        printf("\nYour field is correct!");
        }
    printf("\nDo you want to save your field? Enter Y or N: ");
    do {scanf("%c",&answer);} while ( getchar() != '\n' );
    if((answer == 'Y') || (answer == 'y')){
        printf("\nEnter a filepath to save the file:");
        do {scanf("%s",filepath);} while ( getchar() != '\n' );
        WriteFlag(filepath,field,sizeX,sizeY);
    }
    printf("\nDo you want to save your field as .vtk File?Enter Y or N: ");
    do {scanf("%c",&answer);} while ( getchar() != '\n' );
    if((answer == 'Y') || (answer == 'y')){
        printf("\nEnter a filepath to save the file:");
        do {scanf("%s",filepath);} while ( getchar() != '\n' );
        writeVTKfileFor2DintegerField(filepath,"shows obstacles",field,sizeX,sizeY,delx,dely);
    }
    printf("\nYou're finished!\n");
}

short CorrectnessCheck(short** flag, int sizeX, int sizeY)
{
    for(int i=2;i<sizeX-1;i++){
        for(int j=2;j<sizeY-1;j++){
            if(flag[i][j]==C_B){
                if(flag[i+1][j]==C_F && flag[i-1][j]==C_F){
                    printf("\nError: (%d,%d) has fluid to the east and west.\n",i,j);
                    return 1;
                }
                if(flag[i][j+1]==C_F && flag[i][j-1]==C_F){
                    printf("\nError: (%d,%d) has fluid to the north and south.\n",i,j);
                    return 2;
                }
            }
        }
    }
    return 0;
}

void EditFlag(short** flag, int sizeX, int sizeY, double dx, double dy)
{
    char run='Y';
    char answer='Y';
    int isize;
    int jsize;
    int ipos;
    int jpos;
    char* filepath;
    char obstacle;
    printf("\nYou're editing an existing field. You can remove and add obstacles.");
    while((run == 'Y') || (run == 'y')){
        printf("\nDo you want to save the field as .vtk for visualization? Enter Y or N");
        do {scanf("%c",&answer);} while ( getchar() != '\n' );
        if((answer == 'Y') || (answer == 'y')){
            printf("\nEnter a filepath to save the file:");
            do {scanf("%s",filepath);} while ( getchar() != '\n' );
            writeVTKfileFor2DintegerField(filepath,"shows obstacles",flag,sizeX,sizeY,dx,dy);
        }
        printf("\nDo you want to add or remove an obstacle? Enter A or R: ");
        do {scanf("%c",&answer);} while ( getchar() != '\n' );
        if((answer == 'A') || (answer == 'a')){
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
            printf("\nWhat type of obstacle do you want to add? Press h for help: ");
            do {scanf("%c",&obstacle);} while ( getchar() != '\n' );
            if(obstacle=='h'|| obstacle == 'H'){
                printf("\nYou can choose between the following types of obstacles:\nR for rectangle\nC for Creative Mode");
                printf("\nPick a type:");
                do {scanf("%c",&obstacle);} while ( getchar() != '\n' );
            }

            switch(obstacle)
            {
            case 'R': case 'r':
                printf("\nSet the position for the lower left block of the rectangle: in x-direction: ");
                do {scanf("%d",&ipos);} while ( getchar() != '\n' );
                printf("\nAnd now in y-direction: ");
                do {scanf("%d",&jpos);} while ( getchar() != '\n' );
                if((ipos<1) | (jpos<1) | ((ipos+isize)>sizeX) | ((jpos+jsize) > sizeY)){
                    printf("\nWrong positions or size too big. No rectangle added.");
                }else{
                    AddRectangle(flag,isize,jsize,ipos,jpos,C_B);
                }
                break;

            default:
                printf("\nNo valid input.");
                break;
            }
        }
        if((answer == 'R') || (answer == 'r')){
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
            printf("\nWhat type of obstacle do you want to remove? Press h for help: ");
            do {scanf("%c",&obstacle);} while ( getchar() != '\n' );
            if(obstacle=='h'|| obstacle == 'H'){
                printf("\nYou can choose between the following types:\nR for rectangle\nC for Creative Mode");
                printf("\nPick a type:");
                do {scanf("%c",&obstacle);} while ( getchar() != '\n' );
            }

            switch(obstacle)
            {
            case 'R': case 'r':
                printf("\nSet the position for the lower left block of the rectangle: in x-direction: ");
                do {scanf("%d",&ipos);} while ( getchar() != '\n' );
                printf("\nAnd now in y-direction: ");
                do {scanf("%d",&jpos);} while ( getchar() != '\n' );
                if((ipos<1) | (jpos<1) | ((ipos+isize)>sizeX) | ((jpos+jsize) > sizeY)){
                    printf("\nWrong positions or size too big. No rectangle removed.");
                }else{
                    AddRectangle(flag,isize,jsize,ipos,jpos,C_F);
                }
                break;

            default:
                printf("\nNo valid input.");
                break;
            }
        }
        printf("Do you want to continue editing? Enter Y or N: ");
        do {scanf("%c",&run);} while ( getchar() != '\n' );
    }
    printf("\nDo you want to save your field? Enter Y or N: ");
    do {scanf("%c",&answer);} while ( getchar() != '\n' );
    if((answer == 'Y') || (answer == 'y')){
        printf("\nEnter a filepath to save the file:");
        do {scanf("%s",filepath);} while ( getchar() != '\n' );
        WriteFlag(filepath,flag,sizeX,sizeY);
    }
    printf("\nYou're finished!\n");
}


void EditFlagFromFile(const char *fileName, double dx, double dy)
{
    int sizeX;
    int sizeY;
    short** flag=ReadFlag(fileName,&sizeX,&sizeY);
    EditFlag(flag,sizeX,sizeY,dx,dy);
}

void AddRectangle(short** flag, int isize,int jsize,int ipos,int jpos, short BorF)
{
    for(int i=ipos;i<ipos+isize;i++){
        for(int j=jpos;j<jpos+jsize;j++){
            flag[i-1][j-1]=BorF;
        }
    }
}
