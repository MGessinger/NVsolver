#include "FlagCreator.h"

void WriteFlag (const char* fileName, char** flag, int sizeX, int sizeY)
{
	FILE *datei;
	datei = fopen(fileName,"wb");
	if (datei == NULL)
	{
		printf("Can't open the file.\n");
	}
	if(datei != NULL){
		fwrite(&sizeX,sizeof(int),1,datei);
		fwrite(&sizeY,sizeof(int),1,datei);
		for(int i = 0;i<sizeX;i++){
			fwrite(flag[i],sizeof(short),sizeY,datei);
		}
	}
	fclose (datei);
}

char** ReadFlag (const char* fileName, int* sizeX, int* sizeY)
{
	FILE* datei;
	datei = fopen(fileName,"rb");
	char **field = NULL;
	if(datei == NULL){
		printf("Can't open the file.\n");
	}
	if (datei != NULL){
		fread(sizeX,sizeof(int),1,datei);
		fread(sizeY,sizeof(int),1,datei);
		field = mallocFlag(*sizeX,*sizeY);
		for(int i = 0;i<*sizeX;i++){
			fread(field[i],sizeof(short),*sizeY,datei);
		}
	}
	fclose(datei);
	return field;
}

char** mallocFlag (int sizeX, int sizeY)
{
	char** field = malloc((sizeX)*(sizeof(char*)+2));
	for (int i = 0; i<sizeX;i++){
		field[i]=malloc(sizeY*(sizeof(char)+2));
	}
	return field;
}

int CreateFlag ()
{
	int sizeX;
	int sizeY;
	printf("First you need to create a field. Enter size in x-direction: ");
	do {scanf("%d",&sizeX);} while ( getchar() != '\n' );
	while(sizeX<1){
		printf("sizeX is too small. Pick a greater number: ");
		do {scanf("%d",&sizeX);} while ( getchar() != '\n' );
	}
	printf("Now enter size in y-direction: ");
	do {scanf("%d",&sizeY);} while ( getchar() != '\n' );
	while(sizeX<1){
		printf("sizeY is too small. Pick a greater number: ");
		do {scanf("%d",&sizeY);} while ( getchar() != '\n' );
	}
	char** field = mallocFlag(sizeX,sizeY);
	for(int i = 0; i<sizeX;i++){
		for(int j = 0;j<sizeX;j++){
			field[i][j]=C_F;
		}
	}
	char run='Y';
	char answer;
	int isize;
	int jsize;
	int ipos;
	int jpos;
	short correct = 0;
	char filepath[100];
	double xlength;
	double ylength;
	char obstacle;
	printf("You created a field. Would you like to add an obstacle? Enter Y or N: ");
	do {scanf("%c",&run);} while ( getchar() != '\n' );
	if((run == 'Y') || (run == 'y')){
		printf("The field has size %d*%d.Please remember the rules for obstacles!\n",sizeX,sizeY);
	}
	printf("We need more Information for Visualization. Please enter xlength: ");
	do {scanf("%lf",&xlength);} while ( getchar() != '\n' );
	printf("And now ylength: ");
	do {scanf("%lf",&ylength);} while ( getchar() != '\n' );
	double delx = xlength/(double)sizeX;
	double dely = ylength/(double)sizeY;
	while((run == 'Y') || (run == 'y')){

		printf("What type of obstacle do you want to add? Press h for help: ");
		do {scanf("%c",&obstacle);} while ( getchar() != '\n' );
		if(obstacle=='h'|| obstacle == 'H'){
			printf("You can choose between the following types of obstacles:\nR for rectangle\nC for Creative Mode");
			printf("\nPick a type:");
			do {scanf("%c",&obstacle);} while ( getchar() != '\n' );
		}

		switch(obstacle)
		{
			case 'R': case 'r':
				printf("Set the position for the lower left block of the rectangle: in x-direction: ");
				do {scanf("%d",&ipos);} while ( getchar() != '\n' );
				printf("And now in y-direction: ");
				do {scanf("%d",&jpos);} while ( getchar() != '\n' );
				printf("Which size should the rectangle have in x-direction? Enter a number: ");
				do {scanf("%d",&isize);} while ( getchar() != '\n' );
				while(isize>sizeX){
					printf("Pick a smaller size please. Enter a new number: ");
					do {scanf("%d",&isize);} while ( getchar() != '\n' );
				}
				printf("Which size should the rectangle have in y-direction? Enter a number: ");
				do {scanf("%d",&jsize);} while ( getchar() != '\n' );
				while(jsize>sizeY){
					printf("Pick a smaller size please. Enter a new number: ");
					do {scanf("%d",&jsize);} while ( getchar() != '\n' );
				}
				if((ipos<1) || (jpos<1) || ((ipos+isize-1)>sizeX) || ((jpos+jsize-1) > sizeY)){
					printf("Wrong positions or size too big. No rectangle added.\n");
				}else{
					AddRectangle(field,isize,jsize,ipos,jpos,C_B);
				}
				break;
			case 'C': case 'c':
				printf("Set the position for the lower left block in x-direction: ");
				do {scanf("%d",&ipos);} while ( getchar() != '\n' );
				printf("And now in y-direction: ");
				do {scanf("%d",&jpos);} while ( getchar() != '\n' );
				if((ipos<1) | (jpos<1)){
					printf("Wrong positions. No CreativeMode possible.\n");
				}else{
					CreativeMode(field,sizeX,sizeY,ipos,jpos);
				}
				break;
			default:
				printf("No valid input.\n");
				break;
		}
		printf("Do you want to check your field in Paraview? Enter Y or N:");
		do {scanf("%c",&answer);} while ( getchar() != '\n' );
		if(answer == 'Y' || answer == 'y'){
			printf("Enter a filepath to save the file:");
			do {scanf("%s",filepath);} while ( getchar() != '\n' );
			lattice *grid = malloc(sizeof(struct lattice));
			grid->delx = delx;
			grid->dely = dely;
			grid->imax = sizeX;
			grid->jmax = sizeY;
			writeVTKfileFor2DintegerField(filepath,"shows_obstacles",field,grid);
			free(grid);
		}
		printf("Do you want to add another obstacle? Enter Y or N: ");
		do {scanf("%c",&run);} while ( getchar() != '\n' );
	}
	printf("You have created a field with obstacles. Do you wanna check if it's correct? Enter Y or N: ");
	do {scanf("%c",&answer);} while ( getchar() != '\n' );
	if((answer == 'Y') || (answer == 'y')){
		correct = CorrectnessCheck(field,sizeX,sizeY);
	}
	while(correct > 0){
		printf("Your field is not correct. Check the error(s) above and edit it.\n");
		EditFlag(field,sizeX,sizeY,delx,dely);
		correct = CorrectnessCheck(field,sizeX,sizeY);
	}
	if((correct == 0 && answer == 'Y')||(correct == 0 && answer == 'y')){
		printf("Your field is correct!\n");
	}
	printf("Do you want to save your field? Enter Y or N: ");
	do {scanf("%c",&answer);} while ( getchar() != '\n' );
	if((answer == 'Y') || (answer == 'y')){
		printf("Enter a filepath to save the file:");
		do {scanf("%s",filepath);} while ( getchar() != '\n' );
		WriteFlag(filepath,field,sizeX,sizeY);
	}
	printf("Do you want to save your field as .vtk File?Enter Y or N: ");
	do {scanf("%c",&answer);} while ( getchar() != '\n' );
	if((answer == 'Y') || (answer == 'y')){
		printf("Enter a filepath to save the file:");
		do {scanf("%s",filepath);} while ( getchar() != '\n' );
		lattice *grid = malloc(sizeof(struct lattice));
		grid->delx = delx;
		grid->dely = dely;
		grid->imax = sizeX;
		grid->jmax = sizeY;
		writeVTKfileFor2DintegerField(filepath,"shows_obstacles",field,grid);
		free(grid);
	}
	printf("\nYou're finished!\n");
	return 0;
}

char CorrectnessCheck (char**flag, int sizeX, int sizeY)
{
	for(int i = 2;i<sizeX-1;i++){
		for(int j = 2;j<sizeY-1;j++){
			if(flag[i][j]==C_B){
				if(flag[i+1][j]==C_F && flag[i-1][j]==C_F){
					printf("\nError: (%d,%d) has fluid to the east and west.\n",i+1,j+1);
					return 1;
				}
				if(flag[i][j+1]==C_F && flag[i][j-1]==C_F){
					printf("\nError: (%d,%d) has fluid to the north and south.\n",i+1,j+1);
					return 2;
				}
			}
		}
	}
	return 0;
}

void EditFlag (char **flag, int sizeX, int sizeY, double dx, double dy)
{
	char run='Y';
	char answer='Y';
	int isize;
	int jsize;
	int ipos;
	int jpos;
	char filepath[100];
	char obstacle;
	printf("\nYou're editing an existing field. You can remove and add obstacles.\n");
	while((run == 'Y') || (run == 'y')){
		printf("Do you want to save the field as .vtk for visualization? Enter Y or N: ");
		do {scanf("%c",&answer);} while ( getchar() != '\n' );
		if((answer == 'Y') || (answer == 'y')){
			printf("Enter a filepath to save the file:");
			do {scanf("%s",filepath);} while ( getchar() != '\n' );
			lattice *grid = malloc(sizeof(struct lattice));
			grid->delx = dx;
			grid->dely = dy;
			grid->imax = sizeX;
			grid->jmax = sizeY;
			writeVTKfileFor2DintegerField(filepath,"shows_obstacles",flag,grid);
			free(grid);
		}
		printf("Do you want to add or remove an obstacle? CreativeMode is available in both. Enter A or R: ");
		do {scanf("%c",&answer);} while ( getchar() != '\n' );
		if((answer == 'A') || (answer == 'a')){

			printf("What type of obstacle do you want to add? Press h for help: ");
			do {scanf("%c",&obstacle);} while ( getchar() != '\n' );
			if(obstacle=='h'|| obstacle == 'H'){
				printf("You can choose between the following types of obstacles:\nR for rectangle\nC for Creative Mode");
				printf("\nPick a type:");
				do {scanf("%c",&obstacle);} while ( getchar() != '\n' );
			}

			switch(obstacle)
			{
				case 'R': case 'r':
					printf("Set the position for the lower left block of the rectangle: in x-direction: ");
					do {scanf("%d",&ipos);} while ( getchar() != '\n' );
					printf("And now in y-direction: ");
					do {scanf("%d",&jpos);} while ( getchar() != '\n' );
					printf("Which size should the new rectangle have in x-direction? Enter a number: ");
					do {scanf("%d",&isize);} while ( getchar() != '\n' );
					while(isize>sizeX){
						printf("Pick a smaller size please. Enter a new number: ");
						do {scanf("%d",&isize);} while ( getchar() != '\n' );
					}
					printf("Which size should the new rectangle have in y-direction? Enter a number: ");
					do {scanf("%d",&jsize);} while ( getchar() != '\n' );
					while(jsize>sizeY){
						printf("Pick a smaller size please. Enter a new number: ");
						do {scanf("%d",&jsize);} while ( getchar() != '\n' );
					}
					if((ipos<1) | (jpos<1) | ((ipos+isize-1)>sizeX) | ((jpos+jsize-1) > sizeY)){
						printf("Wrong positions or size too big. No rectangle added.\n");
					}else{
						AddRectangle(flag,isize,jsize,ipos,jpos,C_B);
					}
					break;
				case 'C': case 'c':

					printf("Set the position for the lower left block in x-direction: ");
					do {scanf("%d",&ipos);} while ( getchar() != '\n' );
					printf("And now in y-direction: ");
					do {scanf("%d",&jpos);} while ( getchar() != '\n' );
					if((ipos<1) | (jpos<1)){
						printf("Wrong positions. No CreativeMode possible.\n");
					}else{
						CreativeMode(flag,sizeX,sizeY,ipos,jpos);
					}
					break;
				default:
					printf("No valid input.\n");
					break;
			}
		}
		if((answer == 'R') || (answer == 'r')){

			printf("What type of obstacle do you want to remove? Press h for help: ");
			do {scanf("%c",&obstacle);} while ( getchar() != '\n' );
			if(obstacle=='h'|| obstacle == 'H'){
				printf("You can choose between the following types:\nR for rectangle\nC for Creative Mode");
				printf("\nPick a type:");
				do {scanf("%c",&obstacle);} while ( getchar() != '\n' );
			}

			switch(obstacle)
			{
				case 'R': case 'r':
					printf("Of which size do you want to remove in x-direction? Enter a number: ");
					do {scanf("%d",&isize);} while ( getchar() != '\n' );
					while(isize>sizeX){
						printf("Pick a smaller size please. Enter a new number: ");
						do {scanf("%d",&isize);} while ( getchar() != '\n' );
					}
					printf("Of which size do you want to remove in y-direction? Enter a number: ");
					do {scanf("%d",&jsize);} while ( getchar() != '\n' );
					while(jsize>sizeY){
						printf("Pick a smaller size please. Enter a new number: ");
						do {scanf("%d",&jsize);} while ( getchar() != '\n' );
					}
					printf("Set the position for the lower left block of the rectangle: in x-direction: ");
					do {scanf("%d",&ipos);} while ( getchar() != '\n' );
					printf("And now in y-direction: ");
					do {scanf("%d",&jpos);} while ( getchar() != '\n' );
					if((ipos<1) || (jpos<1) || ((ipos+isize-1)>sizeX) || ((jpos+jsize-1) > sizeY)){
						printf("Wrong positions or size too big. No rectangle removed.\n");
					}else{
						AddRectangle(flag,isize,jsize,ipos,jpos,C_F);
					}
					break;
				case 'C': case 'c':
					printf("Set the position for the lower left block in x-direction: ");
					do {scanf("%d",&ipos);} while ( getchar() != '\n' );
					printf("And now in y-direction: ");
					do {scanf("%d",&jpos);} while ( getchar() != '\n' );
					if((ipos<1) | (jpos<1)){
						printf("Wrong positions. No CreativeMode possible.\n");
					}else{
						CreativeMode(flag,sizeX,sizeY,ipos,jpos);
					}
					break;
				default:
					printf("No valid input.\n");
					break;
			}
		}
		printf("Do you want to continue editing? Enter Y or N: ");
		do {scanf("%c",&run);} while ( getchar() != '\n' );
	}
	printf("\nYou're finished with editing!\n");
}

void EditFlagFromFile ()
{
	char fileName[100];
	double xlength;
	double ylength;
	printf("Which file do you want to edit? Enter the filepath:");
	do {scanf("%s",fileName);} while ( getchar() != '\n' );
	printf("We need more Information for Visualization. Please enter xlength: ");
	do {scanf("%lf",&xlength);} while ( getchar() != '\n' );
	printf("And now ylength: ");
	do {scanf("%lf",&ylength);} while ( getchar() != '\n' );
	int sizeX;
	int sizeY;
	char** flag = ReadFlag(fileName,&sizeX,&sizeY);
	double delx = xlength/(double)sizeX;
	double dely = ylength/(double)sizeY;
	EditFlag(flag,sizeX,sizeY,delx,dely);
}

void AddRectangle (char **flag, int isize, int jsize, int ipos, int jpos, short BorF)
{
	for(int i = ipos;i<ipos+isize;i++){
		for(int j = jpos;j<jpos+jsize;j++){
			flag[i-1][j-1]=BorF;
		}
	}
}

void CreativeMode (char**flag, int sizeX, int sizeY, int ipos, int jpos)
{
	char input[101];
	const char trenn[]=".\n";
	char borf;
	int num;
	int jh = jpos;
	int ih = ipos;
	printf("CreativeMode is very dangerous. Please read readme.txt first. Enter your code now:\n");
	do {scanf("%s",input);} while ( getchar() != '\n' );
	char* todo = strtok(input,trenn);
	while(todo != NULL){
		borf = *todo;
		todo = strtok(NULL, trenn);
		num = atoi(todo);
		todo = strtok(NULL, trenn);
		if(((ih+num)>sizeX)||(jh>sizeY)){
			printf("Input error: Try again");
		}else{
			switch(borf)
			{
				case 'b': case 'B':
					for(int i = ih;i<ih+num;i++){
						flag[i-1][jh-1]=C_B;
					}
					ih = ih+num;
					break;
				case 'f': case 'F':
					for (int i = ih;i<ih+num;i++){
						flag[i-1][jh-1]=C_F;
					}
					ih = ih+num;
					break;
				case 'n': case 'N':
					jh++;
					ih = ipos;
					break;
			}
		}
	}
}
