#include <stdio.h>
#include <stdlib.h>
#include <math.h>
float PI=3.14159265359;

/*construct a structure of BMP header*/
typedef struct Bmpheader{
    unsigned short identifier; // 0x0000
    unsigned int filesize; // 0x0002
    unsigned short reserved; // 0x0006
    unsigned short reserved2;
    unsigned int bitmap_dataoffset; // 0x000A
    unsigned int bitmap_headersize; // 0x000E
    unsigned int width; // 0x0012
    unsigned int height; // 0x0016
    unsigned short planes; // 0x001A
    unsigned short bits_perpixel; // 0x001C
    unsigned int compression; // 0x001E
    unsigned int bitmap_datasize; // 0x0022
    unsigned int hresolution; // 0x0026
    unsigned int vresolution; // 0x002A
    unsigned int usedcolors; // 0x002E
    unsigned int importantcolors; // 0x0032
    unsigned int palette; // 0x0036
} Bitmap;

/*construct a structure of RGB*/
typedef struct RGB{
    int R;
    int G;
    int B;
} ImgRGB;

typedef struct YIQ{
    float Y;
    float I;
    float Q;
} ImgYIQ;
typedef struct DCT_YIQ{
    int Y;
    int I;
    int Q;
} DCT_YIQ;
typedef struct data_zig{
    int num;
    int count;
    int code[30];
    int len;
} data_zig;


void inverse_zig(int *data_Y, int *data_I, int *data_Q, int size_Y, int size_I, int size_Q, Bitmap bmpheader);
Bitmap readheader(FILE* fp);
ImgRGB** malloc_2D(int row, int col);
ImgYIQ** malloc_2D_YIQ(int row, int col);
DCT_YIQ** malloc_2D_DCT_YIQ(int row, int col);
void InputData(FILE* fp,ImgRGB **array,int H,int W);
void output_bmp(ImgRGB **RGB,FILE* outfile,Bitmap bmpheader);
void output_bmp_YIQ(ImgYIQ **YIQ,FILE* outfile,Bitmap bmpheader);
void convert_YIQ(ImgRGB **RGB, ImgYIQ **YIQ, int row, int col);
void inverse_YIQ(ImgRGB **RGB, ImgYIQ **YIQ, int row, int col);
void shift_value(ImgYIQ **YIQ, int row, int col);
void in_shift_value(ImgYIQ **YIQ, int row, int col);
ImgYIQ** DFT(ImgYIQ **data);
ImgYIQ** inverse_DFT(ImgYIQ  **data);
DCT_YIQ** blocking(ImgYIQ **data, Bitmap bmpheader);
DCT_YIQ** quantize(ImgYIQ **data);
ImgYIQ** in_quantize(DCT_YIQ **data);
int DC_dif(int value, int f_value);
int* zigzag(DCT_YIQ **data, int height, int width, Bitmap bmpheader, int *zig_Y, int *zig_I, int *zig_Q);
void huffman_Y(int *data, int size);
void huffman_I(int *data, int size);
void huffman_Q(int *data, int size);
void codebook(data_zig *data,int *code_seq, int size, FILE *fp);
int* code_sequence(data_zig *data, int size);
void encode(int *zig_Y, int *zig_I, int *zig_Q, int size_Y, int size_I, int size_Q,  Bitmap bmpheader);
int read_codebook(data_zig *data, FILE *fp);
int encode_generator(int *zig_data, data_zig *codebook, int codebook_len, int zig_size, int *encode_data);
void padding(int *data);
void decode(int *data_Y, int *data_I, int *data_Q, data_zig *codebook_Y, data_zig *codebook_I, data_zig *codebook_Q, Bitmap bmpheader, int size_Y, int size_I, int size_Q);
void encode_data(int *data_Y, int *data_I, int *data_Q, FILE *fp, int size_Y, int size_I, int size_Q);

int main(int argc,char *argv[]){

    FILE *fp=fopen(argv[1],"rb");
    //FILE *fw=fopen(argv[2],"wb");
    int i, j, k, h;

    Bitmap bmpheader=readheader(fp);
    printf("Read header success\n");
    ImgRGB **Data_RGB=malloc_2D(bmpheader.height, bmpheader.width);
    printf("build data_RGB success\n");
    InputData(fp,Data_RGB,bmpheader.height,bmpheader.width);
    printf("reading data sucess data saved into Dat_RGB\n");
    ImgYIQ **Data_YIQ=malloc_2D_YIQ(bmpheader.height, bmpheader.width);
    printf("Build YIQ success\n");
    ImgRGB **in_data_rgb = malloc_2D(bmpheader.height, bmpheader.width);
    printf("inverse data build\n");
    convert_YIQ(Data_RGB, Data_YIQ, bmpheader.height, bmpheader.width);
    printf("Convert YIQ success\n");
    shift_value(Data_YIQ,bmpheader.height, bmpheader.width);
    in_shift_value(Data_YIQ, bmpheader.height, bmpheader.width);
    printf("Shift value success\n");
    // inverse_YIQ(in_data_rgb, Data_YIQ, bmpheader.height, bmpheader.width);
    // printf("inverse done\n");

    //output_bmp(in_data_rgb,fw,bmpheader);

    ///////////////////////////////////////blocking////////////////////////////////////////////////
    DCT_YIQ** zigzag_data=malloc_2D_DCT_YIQ(bmpheader.height, bmpheader.width);
    zigzag_data = blocking(Data_YIQ, bmpheader);
    //output_bmp_YIQ(Data_YIQ, fw, bmpheader);

    // //////////////////////////////////////recover file/////////////////////////////////////////////
    // inverse_YIQ(in_data_rgb,Data_YIQ, bmpheader.height, bmpheader.width);
    // output_bmp(in_data_rgb,fw,bmpheader);
    // printf("output_bmp\n");
    // fclose(fp);

    // ///////////////////////////////////////zigzag/////////////////////////////////////////////////
    int *zig_Y = calloc(bmpheader.height*bmpheader.width, sizeof(int));
    int *zig_I = calloc(bmpheader.height*bmpheader.width, sizeof(int));
    int *zig_Q = calloc(bmpheader.height*bmpheader.width, sizeof(int));
    int *zigzag_size;
    zigzag_size = zigzag(zigzag_data, bmpheader.height, bmpheader.width, bmpheader, zig_Y, zig_I, zig_Q);
    //printf("\n\n\nY:%d I:%d Q:%d\n\n\n", zigzag_size[0], zigzag_size[1], zigzag_size[2]);
    printf("Zigzag done\n");
    fclose(fp);
    /////////////////////////////////Encode////////////////////////////////////////////
    encode(zig_Y, zig_I, zig_Q, zigzag_size[0], zigzag_size[1], zigzag_size[2], bmpheader);
    for(i=0;i<bmpheader.height*bmpheader.width;i++){
        free(Data_RGB[i]);
        free(Data_YIQ[i]);
        free(in_data_rgb[i]);
        free(zigzag_data[i]);
    }
    //fclose(fw);
    free(Data_RGB);
    free(Data_YIQ);
    free(in_data_rgb);
    free(zigzag_data);
    free(zig_Y);
    free(zig_I);
    free(zig_Q);
    
    return 0;
}


/*read header*/
Bitmap readheader(FILE* fp){
    Bitmap x;
    fread(&x.identifier,sizeof(x.identifier),1,fp);
    fread(&x.filesize,sizeof(x.filesize),1,fp);
    fread(&x.reserved,sizeof(x.reserved),1,fp);
    fread(&x.reserved2,sizeof(x.reserved2),1,fp);
    fread(&x.bitmap_dataoffset,sizeof(x.bitmap_dataoffset),1,fp);
    fread(&x.bitmap_headersize,sizeof(x.bitmap_headersize),1,fp);
    fread(&x.width,sizeof(x.width),1,fp);
    fread(&x.height,sizeof(x.height),1,fp);
    fread(&x.planes,sizeof(x.planes),1,fp);
    fread(&x.bits_perpixel,sizeof(x.bits_perpixel),1,fp);
    fread(&x.compression,sizeof(x.compression),1,fp);
    fread(&x.bitmap_datasize,sizeof(x.bitmap_datasize),1,fp);
    fread(&x.hresolution,sizeof(x.hresolution),1,fp);
    fread(&x.vresolution,sizeof(x.vresolution),1,fp);
    fread(&x.usedcolors,sizeof(x.usedcolors),1,fp);
    fread(&x.importantcolors,sizeof(x.importantcolors),1,fp);
    printf("width : %d   ", x.width);
    printf("height : %d\n", x.height);
    return x;
}

/*input data without header into RGB structure*/
void InputData(FILE* fp,ImgRGB **array,int H, int W){
    int temp, i, j;
    for(i=0;i<H;i++){
        for(j=0;j<W;j++){
            temp=fgetc(fp);
            array[i][j].B=temp;
            //printf("%d", array[i][j].B);
            temp=fgetc(fp);
            array[i][j].G=temp;
            //printf("%d", array[i][j].G);
            temp=fgetc(fp);
            array[i][j].R=temp;
            //printf("%d", array[i][j].R);
            //printf("\n");
        }
    }
}

/* A function of making two dimensions memory locate array*/
ImgRGB** malloc_2D(int row, int col){
    ImgRGB **Array, *Data;
    int i;
    Array=(ImgRGB**)malloc(row*sizeof(ImgRGB *));
    Data=(ImgRGB*)malloc(row*col*sizeof(ImgRGB));
    for(i=0; i<row; i++,Data+=col){
        Array[i] = Data;
    }
    return Array;
}

/*output header and data*/
void output_bmp(ImgRGB **RGB,FILE* outfile,Bitmap bmpheader){
    fwrite(&bmpheader.identifier, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.filesize, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.reserved, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.reserved2, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.bitmap_dataoffset, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.bitmap_headersize, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.width, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.height, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.planes, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.bits_perpixel, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.compression, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.bitmap_datasize, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.hresolution, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.vresolution, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.usedcolors, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.importantcolors, sizeof(int), 1 , outfile);
    int x, y;
    for (x=0; x<bmpheader.height; x++){  
        for(y=0; y<bmpheader.width; y++){   
            fwrite(&RGB[x][y].B, sizeof(char),1,outfile);
            fwrite(&RGB[x][y].G, sizeof(char),1,outfile);
            fwrite(&RGB[x][y].R, sizeof(char),1,outfile);
        }
    }
    //printf("x=%d ", x);
}

void output_bmp_YIQ(ImgYIQ **YIQ,FILE* outfile,Bitmap bmpheader){
    fwrite(&bmpheader.identifier, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.filesize, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.reserved, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.reserved2, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.bitmap_dataoffset, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.bitmap_headersize, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.width, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.height, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.planes, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.bits_perpixel, sizeof(short), 1 , outfile);
    fwrite(&bmpheader.compression, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.bitmap_datasize, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.hresolution, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.vresolution, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.usedcolors, sizeof(int), 1 , outfile);
    fwrite(&bmpheader.importantcolors, sizeof(int), 1 , outfile);
    int x, y;
    for (x=0; x<bmpheader.height; x++){  
        for(y=0; y<bmpheader.width; y++){   
            fwrite(&YIQ[x][y].Y, sizeof(char),1,outfile);
            fwrite(&YIQ[x][y].I, sizeof(char),1,outfile);
            fwrite(&YIQ[x][y].Q, sizeof(char),1,outfile);
        }
    }
}

ImgYIQ** malloc_2D_YIQ(int row, int col){
    //printf("Building YIQ array...\n");
    ImgYIQ **Array, *Data;
    int i;
    Array=(ImgYIQ**)malloc(row*sizeof(ImgYIQ *));
    Data=(ImgYIQ*)malloc(row*col*sizeof(ImgYIQ));
    for(i=0; i<row; i++,Data+=col){
        Array[i] = Data;
    }
    return Array;
}

DCT_YIQ** malloc_2D_DCT_YIQ(int row, int col){
    //printf("Building YIQ array...\n");
    DCT_YIQ **Array, *Data;
    int i;
    Array=(DCT_YIQ**)malloc(row*sizeof(DCT_YIQ *));
    Data=(DCT_YIQ*)malloc(row*col*sizeof(DCT_YIQ));
    for(i=0; i<row; i++,Data+=col){
        Array[i] = Data;
    }
    return Array;
}

void convert_YIQ(ImgRGB **RGB, ImgYIQ **YIQ, int row, int col){
    printf("Convert YIQ...\n");
    int i,j;
    for(i=0;i<row;i++){
        for(j=0;j<col;j++){
            YIQ[i][j].Y = 0.299*RGB[i][j].R + 0.587*RGB[i][j].G + 0.114*RGB[i][j].B;
            YIQ[i][j].I = 0.596*RGB[i][j].R - 0.274*RGB[i][j].G - 0.322*RGB[i][j].B;
            YIQ[i][j].Q = 0.211*RGB[i][j].R - 0.523*RGB[i][j].G + 0.312*RGB[i][j].B;
        }
    }
}

void inverse_YIQ(ImgRGB **RGB, ImgYIQ **YIQ, int row, int col){
    printf("Inverse YIQ...\n");
    int i,j;
    for(i=0;i<row;i++){
        for(j=0;j<col;j++){
            RGB[i][j].R = YIQ[i][j].Y + 0.956*YIQ[i][j].I + 0.619*YIQ[i][j].Q;
            RGB[i][j].G = YIQ[i][j].Y - 0.272*YIQ[i][j].I - 0.647*YIQ[i][j].Q;
            RGB[i][j].B = YIQ[i][j].Y - 1.106*YIQ[i][j].I + 1.703*YIQ[i][j].Q;
        }
    }
}

void shift_value(ImgYIQ **YIQ, int row, int col){
    int x,y;
    for(x=0;x<row;x++){
        for(y=0;y<col;y++){
            YIQ[x][y].Y = YIQ[x][y].Y-128.0;
            YIQ[x][y].I = YIQ[x][y].I-128.0;
            YIQ[x][y].Q = YIQ[x][y].Q-128.0;
        }
    }
}

void in_shift_value(ImgYIQ **YIQ, int row, int col){
    int x,y;
    for(x=0;x<row;x++){
        for(y=0;y<col;y++){
            YIQ[x][y].Y = YIQ[x][y].Y+128.0;
            YIQ[x][y].I = YIQ[x][y].I+128.0;
            YIQ[x][y].Q = YIQ[x][y].Q+128.0;
        }
    }
}

DCT_YIQ** blocking(ImgYIQ **data, Bitmap bmpheader){
    //FILE *DFT_pic=fopen("inverse.bmp", "wb");
    int i,j,k,h,flag;
    int block_h=bmpheader.height/8;
    int block_w=bmpheader.width/8;
    printf("block_h = %d  block_w = %d\n", block_h, block_w);
    int remain_h=bmpheader.height%8;
    int remain_w=bmpheader.width%8;
    printf("remain height = %d  remain width = %d\n", remain_h, remain_w);
    //ImgRGB **in_data_rgb=malloc_2D(bmpheader.height,bmpheader.width);
    ImgYIQ **DFT_data=malloc_2D_YIQ(8,8);
    ImgYIQ **DFT_out=malloc_2D_YIQ(8,8);
    DCT_YIQ **quant=malloc_2D_DCT_YIQ(8,8);
    DCT_YIQ **DFTy=malloc_2D_DCT_YIQ(bmpheader.height, bmpheader.width);
    DCT_YIQ **in_quant=malloc_2D_DCT_YIQ(8,8);
    ImgYIQ **in_DFT=malloc_2D_YIQ(8,8);
    ImgYIQ **DFTx=malloc_2D_YIQ(bmpheader.height, bmpheader.width);
    DCT_YIQ first_DC;
    int pre_Y, pre_I, pre_Q;
    printf("DFT -> Quantilizing -> Store DC value\n");
    for(i=0;i<block_h;i++){
        for(j=0;j<block_w;j++){
            for(k=0;k<8;k++){
                for(h=0;h<8;h++){
                    DFT_data[k][h].Y = data[i*8+k][j*8+h].Y;
                    DFT_data[k][h].I = data[i*8+k][j*8+h].I;
                    DFT_data[k][h].Q = data[i*8+k][j*8+h].Q;
                }
            }
            
            DFT_out = DFT(DFT_data);   
            quant = quantize(DFT_out);   

            if(i==0 && j==0){
                first_DC.Y = quant[0][0].Y;
                printf("first_DC_Y = %d ", first_DC.Y);
                first_DC.I = quant[0][0].I;
                printf("first_DC_I = %d ", first_DC.I);
                first_DC.Q = quant[0][0].Q; 
                printf("first_DC_Q = %d \n", first_DC.Q);
            }else{
                pre_Y = quant[0][0].Y;
                pre_I = quant[0][0].I;
                pre_Q = quant[0][0].Q;
                quant[0][0].Y = DC_dif(quant[0][0].Y, first_DC.Y);
                quant[0][0].I = DC_dif(quant[0][0].I, first_DC.I);
                quant[0][0].Q = DC_dif(quant[0][0].Q, first_DC.Q);
                first_DC.Y = pre_Y;
                first_DC.I = pre_I;
                first_DC.Q = pre_Q;
            }

            
            for(k=0;k<8;k++){
                for(h=0;h<8;h++){
                    DFTy[i*8+k][j*8+h].Y = quant[k][h].Y;
                    DFTy[i*8+k][j*8+h].I = quant[k][h].I;
                    DFTy[i*8+k][j*8+h].Q = quant[k][h].Q;
                }
            }
        }
    }
    return DFTy;
}

ImgYIQ** DFT(ImgYIQ **data){
    //printf("DFT data\n");
    ImgYIQ **pixel = malloc_2D_YIQ(8,8);
    int u, v, r, s;
    float Cu, Cv, x=0.0, y=0.0, z=0.0;
    for(u=0;u<8;u++){
        for(v=0;v<8;v++){
            pixel[u][v].Y = 0.0;
            pixel[u][v].I = 0.0;
            pixel[u][v].Q = 0.0;
        }
    }
    for(u=0;u<8;u++){
        for(v=0;v<8;v++){
            if(u==0){ Cu=0.707106781; }
            else{ Cu=1.0; }
            if(v==0){ Cv=0.707106781; }
            else{ Cv=1.0;}
            //float x = sqrt(u*v);
            for(r=0;r<8;r++){
                for(s=0;s<8;s++){
                    x+=data[r][s].Y*(cos(((2*r)+1)*u*PI/16))*(cos(((2*s)+1)*v*PI/16));
                    y+=data[r][s].I*(cos(((2*r)+1)*u*PI/16))*(cos(((2*s)+1)*v*PI/16));
                    z+=data[r][s].Q*(cos(((2*r)+1)*u*PI/16))*(cos(((2*s)+1)*v*PI/16));
                }
            }
            pixel[u][v].Y = Cu*Cv*x/4;
            pixel[u][v].I = Cu*Cv*y/4;
            pixel[u][v].Q = Cu*Cv*z/4;
            x=0.0;
            y=0.0;
            z=0.0;
        }
    }
    return pixel;
}

ImgYIQ** inverse_DFT(ImgYIQ **data){
    //printf("Inverse DFT...\n");
    ImgYIQ **inverse_dft=malloc_2D_YIQ(8, 8);
    float Cv, Cu, Y=0.0, I=0.0, Q=0.0;
    int r, s, u, v;
    for(u=0;u<8;u++){
        for(v=0;v<8;v++){
            inverse_dft[u][v].Y = 0.0;
            inverse_dft[u][v].I = 0.0;
            inverse_dft[u][v].Q = 0.0;
        }
    }
    for(r=0;r<8;r++){
        for(s=0;s<8;s++){
            for(u=0;u<8;u++){
                for(v=0;v<8;v++){
                    if(u==0){ Cu=0.707106781; }
                    else{ Cu=1.0; }
                    if(v==0){ Cv=0.707106781; }
                    else{ Cv=1.0; }
                    Y += Cu*Cv*data[u][v].Y*(cos(((2*r)+1)*u*PI*16))*(cos(((2*s)+1)*v*PI/16));
                    I += Cu*Cv*data[u][v].I*(cos(((2*r)+1)*u*PI*16))*(cos(((2*s)+1)*v*PI/16));
                    Q += Cu*Cv*data[u][v].Q*(cos(((2*r)+1)*u*PI*16))*(cos(((2*s)+1)*v*PI/16));
                }
            }

            inverse_dft[r][s].Y= Y/4;
            inverse_dft[r][s].I= I/4;
            inverse_dft[r][s].Q= Q/4;
            Y=0.0;
            I=0.0;
            Q=0.0;
        }
    }
    return inverse_dft;
}

DCT_YIQ** quantize(ImgYIQ **data){
    int quantize_table_Y[8][8] = {
        {16, 11, 10, 16, 24, 40, 51, 61},
        {12, 12, 14, 19, 26, 58, 60, 55},
        {14, 13, 16, 24, 40, 57, 69, 56},
        {14, 17, 22, 29, 51, 87, 80, 62},
        {18, 22, 37, 56, 68,109,103, 77},
        {24, 35, 55, 64, 81,104,113, 92},
        {49, 64, 78, 87,103,121,120,101},
        {72, 92, 95, 98,112,100,103, 99}
    };
    int quantize_table_IQ[8][8] = {
        {17,18,24,47,99,99,99,99},
        {18,21,26,66,99,99,99,99},
        {24,26,56,99,99,99,99,99},
        {47,66,99,99,99,99,99,99},
        {99,99,99,99,99,99,99,99},
        {99,99,99,99,99,99,99,99},
        {99,99,99,99,99,99,99,99},
        {99,99,99,99,99,99,99,99}
    };

    DCT_YIQ **quant=malloc_2D_DCT_YIQ(8,8);
    int i,j;
    for(i=0;i<8;i++){
        for(j=0;j<8;j++){
            quant[i][j].Y = data[i][j].Y/quantize_table_Y[i][j];
            //printf("quant: %d  %f\n", quant[i][j].Y, data[i][j].Y);
            quant[i][j].I = data[i][j].I/quantize_table_IQ[i][j];
            quant[i][j].Q = data[i][j].Q/quantize_table_IQ[i][j];
        }
    }
    return quant;
}

ImgYIQ** in_quantize(DCT_YIQ **data){
    float quantize_table_Y[8][8] = {
        {16.0,11.0,10.0,16.0,24.0,40.0,51.0,61.0},
        {12.0,12.0,14.0,19.0,26.0,58.0,60.0,55.0},
        {14.0,13.0,16.0,24.0,40.0,57.0,69.0,56.0},
        {14.0,17.0,22.0,29.0,51.0,87.0,80.0,62.0},
        {18.0,22.0,37.0,56.0,68.0,109.0,103.0,77.0},
        {24.0,35.0,55.0,64.0,81.0,104.0,113.0,92.0},
        {49.0,64.0,78.0,87.0,103.0,121.0,120.0,101.0},
        {72.0,92.0,95.0,98.0,112.0,100.0,103.0,99.0}
    };
    float quantize_table_IQ[8][8] = {
        {17.0,18.0,24.0,47.0,99.0,99.0,99.0,99.0},
        {18.0,21.0,26.0,66.0,99.0,99.0,99.0,99.0},
        {24.0,26.0,56.0,99.0,99.0,99.0,99.0,99.0},
        {47.0,66.0,99.0,99.0,99.0,99.0,99.0,99.0},
        {99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0},
        {99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0},
        {99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0},
        {99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0}
    };

    ImgYIQ **quant=malloc_2D_YIQ(8,8);
    int i,j;
    for(i=0;i<8;i++){
        for(j=0;j<8;j++){
            quant[i][j].Y = data[i][j].Y*quantize_table_Y[i][j];
            //printf("in_quant: %f  %d\n", quant[i][j].Y, data[i][j].Y);
            quant[i][j].I = data[i][j].I*quantize_table_IQ[i][j];
            quant[i][j].Q = data[i][j].Q*quantize_table_IQ[i][j];
        }
    }
    return quant;
}

int DC_dif(int value, int f_value){
    int x;
    x = value-f_value;
    return x;
}

int* zigzag(DCT_YIQ **data, int height, int width, Bitmap bmpheader,int *zig_Y, int *zig_I, int *zig_Q){
    int i, j , k, h, flag_Y=0, flag_I=0, flag_Q=0, Y=0, I=0, Q=0;
    DCT_YIQ **num = malloc_2D_DCT_YIQ(8,8);
    int *size = calloc(3, sizeof(int));
    int *zigzag_Y = calloc(height*width, sizeof(int));
    int *zigzag_I = calloc(height*width, sizeof(int));
    int *zigzag_Q = calloc(height*width, sizeof(int));
    int x[]={0,0,1,2,1,0,0,1,2,3,4,3,2,1,0,0,1,2,3,4,5,6,5,4,3,2,1,0,0,1,2,3,4,5,6,7,7,6,5,
            4,3,2,1,2,3,4,5,6,7,7,6,5,4,3,4,5,6,7,7,6,5,6,7,7};
    int y[]={0,1,0,0,1,2,3,2,1,0,0,1,2,3,4,5,4,3,2,1,0,0,1,2,3,4,5,6,7,6,5,4,3,2,1,0,1,2,3,
            4,5,6,7,7,6,5,4,3,2,3,4,5,6,7,7,6,5,4,5,6,7,7,6,7};
    for(i=0;i<height/8;i++){
        for(j=0;j<width/8;j++){
            for(k=0;k<8;k++){
                for(h=0;h<8;h++){
                    num[k][h].Y = data[i*8+k][j*8+h].Y;
                    num[k][h].I = data[i*8+k][j*8+h].I;
                    num[k][h].Q = data[i*8+k][j*8+h].Q;
                }
            }
            
            for(k=0;k<64;k++){
                if(k==63){
                    *(zigzag_Y+Y) = flag_Y;
                    *(zigzag_Y+Y+1) = num[x[k]][y[k]].Y;
                    Y=Y+2;
                    flag_Y=0;
                    *(zigzag_I+I) = flag_I;
                    *(zigzag_I+I+1) = num[x[k]][y[k]].I;
                    I=I+2;
                    flag_I=0;
                    *(zigzag_Q+Q) = flag_Q;
                    *(zigzag_Q+Q+1) = num[x[k]][y[k]].Q;
                    Q=Q+2;
                    flag_Q=0;
                    break;
                }
                if(num[x[k]][y[k]].Y != 0){
                    *(zigzag_Y+Y) = flag_Y;
                    *(zigzag_Y+Y+1) = num[x[k]][y[k]].Y;
                    Y=Y+2;
                    flag_Y=0;
                }else{
                    flag_Y++;
                }
                if(num[x[k]][y[k]].I != 0){
                    *(zigzag_I+I) = flag_I;
                    *(zigzag_I+I+1) = num[x[k]][y[k]].I;
                    I=I+2;
                    flag_I=0;
                }else{
                    flag_I++;
                }
                if(num[x[k]][y[k]].Q != 0){
                    *(zigzag_Q+Q) = flag_Q;
                    *(zigzag_Q+Q+1) = num[x[k]][y[k]].Q;
                    Q=Q+2;
                    flag_Q=0;
                }else{
                    flag_Q++;
                }
                
            }
        }
    }
    //inverse_zig(zigzag_Y, zigzag_I, zigzag_Q, Y, I, Q, bmpheader);
    for(i=0;i<Y;i++){
        zig_Y[i] = zigzag_Y[i];
    }
    size[0]=Y;
    for(i=0;i<I;i++){
        zig_I[i] = zigzag_I[i];
    }
    size[1]=I;
    for(i=0;i<Q;i++){
        zig_Q[i] = zigzag_Q[i];
    }
    size[2]=Q;
    huffman_Y(zigzag_Y, Y);
    huffman_I(zigzag_I, I);
    huffman_Q(zigzag_Q, Q);
    return size;
}

void huffman_Y(int *data, int size){
    printf("\nHuffman Y processing...\n");
    printf("size of  zigzag Y %d\n", size);
    FILE *fp = fopen("codebook_Y.txt", "wb");
    data_zig *zig_data = calloc(3000, sizeof(data_zig));
    int *code_seq = calloc(1000, sizeof(int));
    int i, j, flag=0, append=0;
    printf("\n");
    for(i=0;i<3000;i++){   //initialize
        zig_data[i].num = 0;
        zig_data[i].count = 0;
    }
    zig_data[flag].num = data[0];
    zig_data[flag].count = 1;
    flag++;
    for(i=1;i<size;i++){
        for(j=0;j<flag;j++){
            if(data[i] == zig_data[j].num){
                zig_data[j].count++;
                append=1;
                break;
            }
        }
        if(append!=1){
            zig_data[flag].num = data[i];
            zig_data[flag].count = 1;
            flag++;
            append=0;
        }else{
            append=0;
        }
    }
    //printf("code sequence Y:\n");
    code_seq = code_sequence(zig_data, flag);
    // for(i=0;i<100;i++){
    //     printf("%d", code_seq[i]);
    // }
    printf("\n");
    int min, wanted;
    int count[flag];
    for(i=0;i<flag;i++){
        count[i]=0;
    }
    data_zig *sort=calloc(flag, sizeof(data_zig));
    for(i=0;i<flag;i++){              //Sorting method             
        min = 1000000;                  
        for(j=0;j<flag;j++){
            if(zig_data[j].count<min && count[j]==0){
                min = zig_data[j].count;
                wanted = j;
            }
        }
        count[wanted]=1;                
        sort[i] = zig_data[wanted];
    }
    // for(i=0;i<flag;i++){
    //     fprintf(fp, "%d   %d\n", sort[i].num, sort[i].count);    
    // }
    codebook(sort, code_seq, flag, fp);
    fclose(fp);
}

void huffman_I(int *data, int size){
    printf("Huffman I processing...\n");
    printf("size of  zigzag I %d\n", size);
    FILE *fp = fopen("codebook_I.txt", "wb");
    data_zig *zig_data = calloc(3000, sizeof(data_zig));
    int *code_seq = calloc(1000, sizeof(int));
    int i, j, flag=0, append=0;
    for(i=0;i<3000;i++){   //initialize
        zig_data[i].num = 0;
        zig_data[i].count = 0;
    }
    zig_data[flag].num = data[0];
    zig_data[flag].count = 1;
    for(i=1;i<size;i++){
        for(j=0;j<flag;j++){
            if(data[i] == zig_data[j].num){
                zig_data[j].count++;
                append=1;
                break;
            }
        }
        if(append!=1){
            //printf("appended_I\n");
            zig_data[flag].num = data[i];
            zig_data[flag].count = 1;
            flag++;
            append=0;
        }else{
            append=0;
        }
    }
    //printf("Code Sequence I:\n");
    code_seq = code_sequence(zig_data, flag);
    // printf("\n");
    int min, wanted;
    int count[flag];
    for(i=0;i<flag;i++){
        count[i]=0;
    }
    data_zig *sort=calloc(flag, sizeof(data_zig));
    for(i=0;i<flag;i++){              //Sorting method             
        min = 1000000;                  
        for(j=0;j<flag;j++){
            if(zig_data[j].count<min && count[j]==0){
                min = zig_data[j].count;
                wanted = j;
            }
        }
        count[wanted]=1;                
        sort[i] = zig_data[wanted];
    }

    // for(i=0;i<flag;i++){
    //     fprintf(fp, "%d   %d\n", sort[i].num, sort[i].count);    
    // }
    codebook(sort, code_seq, flag, fp);
    fclose(fp);
}

void huffman_Q(int *data, int size){
    printf("Huffman Q processing...\n");
    printf("size of zigzag Q %d\n", size);
    FILE *fp = fopen("codebook_Q.txt", "wb");
    data_zig *zig_data = calloc(3000, sizeof(data_zig));
    int *code_seq = calloc(1000, sizeof(int));
    int i, j, flag=0, append=0;
    for(i=0;i<3000;i++){   //initialize
        zig_data[i].num = 0;
        zig_data[i].count = 0;
    }
    zig_data[flag].num = data[0];
    zig_data[flag].count = 1;
    for(i=1;i<size;i++){
        for(j=0;j<flag;j++){
            if(data[i] == zig_data[j].num){
                zig_data[j].count++;
                append=1;
                break;
            }
        }
        if(append!=1){
            zig_data[flag].num = data[i];
            zig_data[flag].count = 1;
            flag++;
            append=0;
        }else{
            append=0;
        }
    }
    //printf("Code Sequence Q:\n");
    code_seq = code_sequence(zig_data, flag);
    //printf("\n");
    int min, wanted;
    int count[flag];
    for(i=0;i<flag;i++){
        count[i]=0;
    }
    data_zig *sort=calloc(flag, sizeof(data_zig));
    for(i=0;i<flag;i++){              //Sorting method             
        min = 1000000;                  
        for(j=0;j<flag;j++){
            if(zig_data[j].count<min && count[j]==0){
                min = zig_data[j].count;
                wanted = j;
            }
        }
        count[wanted]=1;                
        sort[i] = zig_data[wanted];
    }
    // for(i=0;i<flag;i++){
    //     fprintf(fp, "%d   %d\n", sort[i].num, sort[i].count);    
    // }
    codebook(sort, code_seq, flag, fp);
    fclose(fp);
}

int* code_sequence(data_zig *data, int size){   //correct
    printf("code sequence generating...\n");
    //FILE *fp=fopen("sort.txt", "wb");
    int i, j, k, min, target, wanted=0, total=0, part=0;
    int sort_count[size];
    int *place = calloc(size-1, sizeof(int));
    int *arr = calloc(size, sizeof(int));
    data_zig sort[size];
    for(i=0;i<size;i++){
        sort_count[i] = 0;
    }
    for(i=0;i<size;i++){     //sort min -> max
        min=1000000;
        for(j=0;j<size;j++){
            if(min>data[j].count && sort_count[j]==0){
                min = data[j].count;
                target = j;
            }    
        }
        total = total + data[target].count;
        sort[i] = data[target];
        sort_count[target] = 1;
        //fprintf(fp, "%d   %d", sort[i].num, sort[i].count);
        arr[i] = data[target].count;
    }

    printf("code sequence length = %d\n", total);
    int kind = size;
    for(i=0;i<size;i++){      //find the 編碼續
        kind = kind-1;            //the size of array will -1 when we add arr[0], arr[1]
        int *arr1 = calloc(kind, sizeof(int));  
        *(arr1) = *(arr) + *(arr+1);    //the new array to save the whole array
        int num = *(arr) + *(arr+1);    
        //printf("%d ", *arr1);
        for(j=1;j<kind;j++){            //assign the whole arr from 2 to the end
            *(arr1+j) = *(arr+j+1);
            //printf("%d ", *(arr1+j));
        }
        //printf("\n");
        int *arr2 = calloc(kind, sizeof(int));
        int *count = calloc(kind, sizeof(int));
        wanted=0;
        for(j=0;j<kind;j++){        //initial count before sorting
            *(count+j)=0;
        }
        for(j=0;j<kind;j++){        // sorting after adding  small -> big
            int min=1000000;
            for(k=0;k<kind;k++){
                if(*(arr1+k)<min && *(count+k)==0){
                    min = *(arr1+k);
                    wanted = k;
                }
                
            }
            *(count+wanted)=1;                
            *(arr2+j) = *(arr1+wanted);
            //printf("%d ", *(arr2+j));
        }
        
        for(j=0;j<kind;j++){           //save the position of arr[0]+arr(1)
            if(*(arr2+j)==num){
                place[part]=j;
                printf("%d",j);
                break;
            }
        }
        part++; 
        arr = arr2;
    }
    //fclose(fp);
    return place;
}

void codebook(data_zig *data, int *code_seq, int size, FILE *fp){
    printf("\nEncoding...\n");
    printf("encode size:%d\n", size);
    int i,j,k;
    int flag=0, x=0;
    for(i=0;i<size;i++){
        data[i].len=0;
        for(j=0;j<25;j++){
            data[i].code[j]=0;
        }
    }
    for(i=0;i<size;i++){  
        //printf("%d  %d ", data[i].num, data[i].count);  
        flag = 0;
        x = i;
        if(x==1 || x==0){   //if it is at position 0 or 1, write code into sort[i].code
            data[i].code[flag] = x;
            //printf("%d", x);
            flag++;
        }
        for(j=0;j<size-1;j++){ //compare the data with 編碼續
            if(x==1 || x==0){
                x = code_seq[j];         
            }else if(x-1<=code_seq[j]){  
                x=x-2;
            }else{
                x=x-1;
            }
            if(x<2){
                data[i].code[flag] = x;  //write code whenever it is at position 1 or 2
                //printf("%d", x);
                flag++;
            }
        }
        //printf("\n");
        data[i].len = flag;
    }
    printf("Done encoding\n");
    for(i=0;i<size;i++){
        fprintf(fp, "%d  ", data[i].num);
        for(j=data[i].len-2;j>=0;j--){
            fprintf(fp, "%d", data[i].code[j]);
        }
        fprintf(fp, "\n");
    }
    printf("Compress txt generated\n");
}

void inverse_zig(int *data_Y, int *data_I, int *data_Q, int size_Y, int size_I, int size_Q, Bitmap bmpheader){
    printf("Into Inverse function\n");
    FILE *fp = fopen("inverse.bmp", "wb");
    printf("bmpheader.width:%d   bmp.height:%d\n", bmpheader.width, bmpheader.height);
    int *in_zigzag_Y = calloc(bmpheader.height*bmpheader.width, sizeof(int));
    int *in_zigzag_I = calloc(bmpheader.height*bmpheader.width, sizeof(int));
    int *in_zigzag_Q = calloc(bmpheader.height*bmpheader.width, sizeof(int));
    DCT_YIQ **zigzag_rebuild = malloc_2D_DCT_YIQ(bmpheader.height,bmpheader.width);
    int position_x[]={0,0,1,2,1,0,0,1,2,3,4,3,2,1,0,0,1,2,3,4,5,6,5,4,3,2,1,0,0,1,2,3,4,5,6,7,7,6,5,
                        4,3,2,1,2,3,4,5,6,7,7,6,5,4,3,4,5,6,7,7,6,5,6,7,7};
    int position_y[]={0,1,0,0,1,2,3,2,1,0,0,1,2,3,4,5,4,3,2,1,0,0,1,2,3,4,5,6,7,6,5,4,3,2,1,0,1,2,3,
                        4,5,6,7,7,6,5,4,3,2,3,4,5,6,7,7,6,5,4,5,6,7,7,6,7};
    int height = bmpheader.height/8;
    int width = bmpheader.width/8;
    int i, j, k, h, flag_Y=0, flag_I=0, flag_Q=0, flag=0;
    int x, y;
    DCT_YIQ **in_quant=malloc_2D_DCT_YIQ(8,8);
    ImgYIQ **in_DFT=malloc_2D_YIQ(8,8);
    ImgYIQ **DFTx=malloc_2D_YIQ(bmpheader.height, bmpheader.width);
    for(i=0;i<bmpheader.height*bmpheader.width;i++){
        in_zigzag_Y[i]=0;
        in_zigzag_I[i]=0;
        in_zigzag_Q[i]=0;
    }
    
    for(i=0;i<size_Y;i=i+2){
        x = data_Y[i];
        y = data_Y[i+1];
        if(x==0){
            in_zigzag_Y[flag_Y] = y;
            flag_Y++;
        }else{
            flag_Y = flag_Y+x;
            in_zigzag_Y[flag_Y] = y;
            flag_Y++;
        }
    }
    
    for(i=0;i<size_I;i=i+2){
        x = data_I[i];
        y = data_I[i+1];
        if(x==0){
            in_zigzag_I[flag_I] = y;
            flag_I++;
        }else{
            flag_I = flag_I+x;
            in_zigzag_I[flag] = y;
            flag_I++;
        }
    }
    //printf("in_zigzag_I done  flag=%d\n", flag_I);
    for(i=0;i<size_Q;i=i+2){
        x = data_Q[i];
        y = data_Q[i+1];
        if(x==0){
            in_zigzag_Q[flag_Q] = y;
            flag_Q++;
        }else{
            flag_Q = flag_Q+x;
            in_zigzag_Q[flag_Q] = y;
            flag_Q++;
        }
    }
    //printf("in_zigzag_Q done flag=%d\n", flag_Q);
    for(i=0;i<height;i++){
        for(j=0;j<width;j++){
            for(k=0;k<64;k++){
                zigzag_rebuild[i*8+position_x[k]][j*8+position_y[k]].Y=in_zigzag_Y[flag];
                zigzag_rebuild[i*8+position_x[k]][j*8+position_y[k]].I=in_zigzag_I[flag];
                zigzag_rebuild[i*8+position_x[k]][j*8+position_y[k]].Q=in_zigzag_Q[flag];
                flag++;
            }
        }
    }
    flag_Y=0;flag_I=0;flag_Q=0;

    for(i=0;i<height;i++){
        for(j=0;j<width;j++){
            for(k=0;k<8;k++){
                for(h=0;h<8;h++){
                    in_quant[k][h].Y = zigzag_rebuild[i*8+k][j*8+h].Y;
                    in_quant[k][h].I = zigzag_rebuild[i*8+k][j*8+h].I;
                    in_quant[k][h].Q = zigzag_rebuild[i*8+k][j*8+h].Q;
                }
            }
            if(i==0 && j==0){
                flag_Y=in_quant[0][0].Y; 
                flag_I=in_quant[0][0].I;
                flag_Q=in_quant[0][0].Q;
            }else{
                in_quant[0][0].Y = in_quant[0][0].Y+flag_Y;
                in_quant[0][0].I = in_quant[0][0].I+flag_I;
                in_quant[0][0].Q = in_quant[0][0].Q+flag_Q;
                flag_Y = in_quant[0][0].Y;
                flag_I = in_quant[0][0].I;
                flag_Q = in_quant[0][0].Q;
            }
            in_DFT = in_quantize(in_quant);
            in_DFT = inverse_DFT(in_DFT);
            for(k=0;k<8;k++){            //initialize
                for(h=0;h<8;h++){
                    DFTx[i*8+k][j*8+h].Y = 0.0;
                    DFTx[i*8+k][j*8+h].I = 0.0;
                    DFTx[i*8+k][j*8+h].Q = 0.0;
                }
            }
            // //printf("data after inverse = %f\n", in_DFT[0][0].Y);
            for(k=0;k<8;k++){
                for(h=0;h<8;h++){
                    DFTx[i*8+k][j*8+h].Y = in_DFT[k][h].Y;
                    DFTx[i*8+k][j*8+h].I = in_DFT[k][h].I;
                    DFTx[i*8+k][j*8+h].Q = in_DFT[k][h].Q;
                }
            }
        }
    }
    ImgRGB **in_data_rgb = malloc_2D(bmpheader.height, bmpheader.width);
    inverse_YIQ(in_data_rgb, DFTx, bmpheader.height, bmpheader.width);
    for(i=0;i<bmpheader.height;i++){
        for(j=0;j<bmpheader.width;j++){
            if(in_data_rgb[i][j].R>255 || in_data_rgb[i][j].G>255 || in_data_rgb[i][j].B>255){
                in_data_rgb[i][j].R=255;
                in_data_rgb[i][j].G=255;
                in_data_rgb[i][j].B=255;
            }
        }
    }
    output_bmp(in_data_rgb, fp,bmpheader);

    free(in_zigzag_Y);
    free(in_zigzag_Q);
    free(in_zigzag_I);
    fclose(fp);    
}

void encode(int *zig_Y, int *zig_I, int *zig_Q, int size_Y, int size_I, int size_Q, Bitmap bmpheader){
    FILE *fp_book_Y=fopen("codebook_Y.txt", "rb");
    FILE *fp_book_I=fopen("codebook_I.txt", "rb");
    FILE *fp_book_Q=fopen("codebook_Q.txt", "rb");
    FILE *fp_out=fopen("compress.txt", "wb");
    int *data_Y = calloc(8000000, sizeof(int));
    int *data_I = calloc(8000000, sizeof(int));
    int *data_Q = calloc(8000000, sizeof(int));
    char trash;
    int i, j, k, h;
    printf("%d %d %d", zig_Y[0], zig_Y[1], zig_Y[2]);
    printf("%d %d %d", zig_I[0], zig_I[1], zig_I[2]);
    printf("%d %d %d", zig_Q[0], zig_Q[1], zig_Q[2]);
    int flag=0, codebooklen_Y, codebooklen_I, codebooklen_Q;
    data_zig *codebook_Y = calloc(200, sizeof(data_zig));
    data_zig *codebook_I = calloc(200, sizeof(data_zig));
    data_zig *codebook_Q = calloc(200, sizeof(data_zig));
    codebooklen_Y = read_codebook(codebook_Y, fp_book_Y);
    codebooklen_I = read_codebook(codebook_I, fp_book_I);
    codebooklen_Q = read_codebook(codebook_Q, fp_book_Q);
    int encode_size_Y, encode_size_I, encode_size_Q;
    encode_size_Y = encode_generator(zig_Y, codebook_Y, codebooklen_Y, size_Y, data_Y);
    printf("encode size:%d\n", encode_size_Y);
    printf("Y encode success\n");
    encode_size_I = encode_generator(zig_I, codebook_I, codebooklen_I, size_I, data_I);
    printf("encode size:%d\n", encode_size_I);
    printf("I encode success\n");
    encode_size_Q = encode_generator(zig_Q, codebook_Q, codebooklen_Q, size_Q, data_Q);
    printf("encode size:%d\n", encode_size_Q);
    printf("Q encode success\n");

    encode_data(data_Y, data_I, data_Q, fp_out,encode_size_Y,encode_size_I,encode_size_Q);
    decode(data_Y, data_I, data_Q, codebook_Y, codebook_I, codebook_Q, bmpheader,encode_size_Y,encode_size_I,encode_size_Q);
    printf("decode done");

    free(data_Y);
    free(data_I);
    free(data_Q);
    free(codebook_Y);
    free(codebook_I);
    free(codebook_Q);
    fclose(fp_book_Y);
    fclose(fp_book_I);
    fclose(fp_book_Q);
    fclose(fp_out);
}

int read_codebook(data_zig *data, FILE *fp){    //corrrect
    printf("reading codebook...\n");
    char trash, x, y;
    int num;
    int flag=0,i;
    while(1){
        x=48;
        y=48;
        fread(&trash, sizeof(char), 1, fp);
        if(feof(fp)) break;
        if(trash == 45){    //numbers < 0
            fread(&x, sizeof(char), 1, fp);
            fread(&y, sizeof(char), 1, fp);
            if(y==32){
                //printf("%c \n", x);
                num=(x-48)*(-1);
                data[flag].num = num;
                fread(&trash, sizeof(char), 1, fp);
            }else{
                //printf("%c %c \n", x, y);
                num=((x-48)*10 + (y-48))*(-1);
                data[flag].num = num;
                fread(&trash, sizeof(char), 1, fp);
                fread(&trash, sizeof(char), 1, fp);
            }
        }else{    //numbers > 0
            x=trash;
            fread(&y, sizeof(char), 1, fp);
            if(y==32){
                //printf("%c \n", x);
                num=(x-48);
                data[flag].num = num;
                fread(&trash, sizeof(char), 1, fp);
            }else{
                //printf("%c %c \n", x, y);
                num=(x-48)*10 + (y-48);
                data[flag].num = num;
                fread(&trash, sizeof(char), 1, fp);
                fread(&trash, sizeof(char), 1, fp);
            }
        }
        //printf("num = %d  ", data[flag].num);
        
        
        for(i=0;i<30;i++){
            fread(&trash, sizeof(char), 1, fp);
            if(trash!=48 && trash!=49){
                break;
            }
            data[flag].code[i]=trash-48;
            //printf("%d", data[flag].code[i]);
        }
        data[flag].len = i;
        //printf("len = %d\n", data[flag].len);
        flag++;
        
    }
    printf("codebook length = %d\n", flag);
    return flag;
}

int encode_generator(int *zig_data, data_zig *codebook, int codebook_len, int zig_size, int *encode_data){
    printf("\nEncode generating\n");
    int *encode = calloc(20000000, sizeof(int));
    int i, j, k, flag=0, app=0;
    //printf("sizeof zigdata: %d\n", zig_size);
    for(i=0;i<zig_size;i++){
        for(j=0;j<codebook_len;j++){
            if(codebook[j].num==zig_data[i]){
                for(k=0;k<codebook[j].len;k++){
                    encode_data[flag] = codebook[j].code[k];
                    flag++;
                }
                break;
            }
        }
    }
    //printf("code length:%d\n", &encode_size);
    return flag;
}

void decode(int *data_Y, int *data_I, int *data_Q, data_zig *codebook_Y, data_zig *codebook_I, data_zig *codebook_Q, Bitmap bmpheader, int size_Y, int size_I, int size_Q){
    printf("\nInverse...\n");
    int i,j,k,h,flag=0, x=0, y=0, z=0;
    int condition=0, state=0;
    int temp[30];
    int *zig_Y = calloc(bmpheader.height*bmpheader.width, sizeof(int));
    int *zig_I = calloc(bmpheader.height*bmpheader.width, sizeof(int));
    int *zig_Q = calloc(bmpheader.height*bmpheader.width, sizeof(int));
    ////////////////////////////Y////////////////////////////////
    for(i=0;i<30;i++){
        temp[i]=0;
    }
    for(i=0;i<size_Y;i++){    
        //printf("Y%d ",i);
        temp[flag] = data_Y[i];
        flag++;
        state=0;
        if(flag==codebook_Y[0].len){
            for(j=0;j<152;j++){
                for(k=0;k<codebook_Y[j].len;k++){
                    if(temp[k] == codebook_Y[j].code[k]){
                        state++;
                    }
                }
                if(state == codebook_Y[j].len){
                    zig_Y[x] = codebook_Y[j].num;
                    x++;
                    condition=1;
                    state=0;
                    break;
                }else{
                    state=0;
                }
            }
            if(condition==1){
                condition=0;
                state=0;
            }
            flag = flag-codebook_Y[j].len;
            for(k=0;k<flag;k++){
                temp[k] = temp[k+codebook_Y[j].len];
            }
        }
    }
    for(j=0;j<152;j++){
        for(k=0;k<codebook_Y[j].len;k++){
            if(temp[k] == codebook_Y[j].code[k]){
                state++;
            }
        }
        if(state == codebook_Y[j].len){
            zig_Y[x] = codebook_Y[j].num;
            x++;
            condition=1;
            state=0;
            break;
        }else{
            state=0;
        }
    }
    //printf("%d  %d  %d  %d  ->%d\n", zig_Y[0], zig_Y[1], zig_Y[2], zig_Y[3], x);
    ///////////////////////////////I///////////////////////////////////////
    flag=0;
    condition=0;
    state=0;
    for(i=0;i<30;i++){
        temp[i]=0;
    }
    for(i=0;i<size_I;i++){   
        //printf("I%d ",i); 
        temp[flag] = data_I[i];
        flag++;
        state=0;
        if(flag==codebook_I[0].len){
            for(j=0;j<57;j++){
                for(k=0;k<codebook_I[j].len;k++){
                    if(temp[k] == codebook_I[j].code[k]){
                        state++;
                    }
                }
                if(state == codebook_I[j].len){
                    zig_I[y] = codebook_I[j].num;
                    y++;
                    condition=1;
                    state=0;
                    break;
                }else{
                    state=0;
                }
            }
            if(condition==1){
                condition=0;
                state=0;
            }
            flag = flag-codebook_I[j].len;
            for(k=0;k<flag;k++){
                temp[k] = temp[k+codebook_I[j].len];
            }
        }
    }
    for(j=0;j<57;j++){
        for(k=0;k<codebook_I[j].len;k++){
            if(temp[k] == codebook_I[j].code[k]){
                state++;
            }
        }
        if(state == codebook_I[j].len){
            zig_I[y] = codebook_I[j].num;
            y++;
            condition=1;
            state=0;
            break;
        }else{
            state=0;
        }
    }
    //printf("%d  %d  %d  %d  ->%d\n", zig_I[0], zig_I[1], zig_I[2], zig_I[3], y);
    /////////////////////////////////Q////////////////////////////////////////
    flag=0;
    condition=0;
    state=0;
    for(i=0;i<30;i++){
        temp[i]=0;
    }
    for(i=0;i<size_Q;i++){   
        //printf("Q%d ",i); 
        temp[flag] = data_Q[i];
        flag++;
        state=0;
        if(flag==codebook_Q[0].len){
            for(j=0;j<35;j++){
                for(k=0;k<codebook_Q[j].len;k++){
                    if(temp[k] == codebook_Q[j].code[k]){
                        state++;
                    }
                }
                if(state == codebook_Q[j].len){
                    zig_Q[z] = codebook_Q[j].num;
                    z++;
                    condition=1;
                    state=0;
                    break;
                }else{
                    state=0;
                }
            }
            if(condition==1){
                condition=0;
                state=0;
            }
            flag = flag-codebook_Q[j].len;
            for(k=0;k<flag;k++){
                temp[k] = temp[k+codebook_Q[j].len];
            }
        }
    }
    for(j=0;j<35;j++){
        for(k=0;k<codebook_Q[j].len;k++){
            if(temp[k] == codebook_Q[j].code[k]){
                state++;
            }
        }
        if(state == codebook_Q[j].len){
            zig_Q[z] = codebook_Q[j].num;
            z++;
            condition=1;
            state=0;
            break;
        }else{
            state=0;
        }
    }
    //printf("%d  %d  %d  %d\n", zig_Q[0], zig_Q[1], zig_Q[2], zig_Q[3]);
    inverse_zig(zig_Y, zig_I, zig_Q, x, y, z, bmpheader);
    printf("\n\ninverse Done\n");
    free(zig_Y);
    free(zig_I);
    free(zig_Q);
}

void encode_data(int *data_Y, int *data_I, int *data_Q, FILE *fp, int size_Y, int size_I, int size_Q){
    printf("Writing Data\n");
    int i,j,k,total;
    int temp[8];
    for(i=0;i<size_Y/8;i++){
        for(j=0;j<8;j++){
            temp[j] = data_Y[i*8+j];
        }
        total =(temp[0])*128 + (temp[1])*64 + (temp[2])*32 + (temp[3])*16 + (temp[4])*8 + (temp[5])*4 + (temp[6])*2 + (temp[7]);
        fprintf(fp, "%c", total);
    }
    for(i=0;i<size_I/8;i++){
        for(j=0;j<8;j++){
            temp[j] = data_I[i*8+j];
        }
        total =(temp[0])*128 + (temp[1])*64 + (temp[2])*32 + (temp[3])*16 + (temp[4])*8 + (temp[5])*4 + (temp[6])*2 + (temp[7]);
        fprintf(fp, "%c", total);
    }
    for(i=0;i<size_Q/8;i++){
        for(j=0;j<8;j++){
            temp[j] = data_Q[i*8+j];
        }
        total =(temp[0])*128 + (temp[1])*64 + (temp[2])*32 + (temp[3])*16 + (temp[4])*8 + (temp[5])*4 + (temp[6])*2 + (temp[7]);
        fprintf(fp, "%c", total);
    }
}