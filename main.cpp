#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mysql.h"

#define LEN 20
#define GENE_LEN 8000000
#define PAM_LEN 20

int check_pam(const char *str,const char *pam){
    for(;*pam;pam++,str++){
        if(*pam=='R' && (*str=='A' || *str=='G')) continue;
        if(*pam=='M' && (*str=='A' || *str=='C')) continue;
        if(*pam=='W' && (*str=='A' || *str=='T')) continue;
        if(*pam=='S' && (*str=='C' || *str=='G')) continue;
        if(*pam=='K' && (*str=='G' || *str=='T')) continue;
        if(*pam=='Y' && (*str=='C' || *str=='T')) continue;
        if(*pam=='H' && (*str=='A' || *str=='C' || *str=='T')) continue;
        if(*pam=='V' && (*str=='A' || *str=='C' || *str=='G')) continue;
        if(*pam=='B' && (*str=='C' || *str=='G' || *str=='T')) continue;
        if(*pam=='D' && (*str=='A' || *str=='G' || *str=='T')) continue;
        if(*pam=='N' || *pam=='n') continue;
        if(*pam==*str) continue;
        return 0;
    }
    return 1;
}

char dna_rev_char(char ch){
    if(ch=='A' || ch=='a') return 'T';
    if(ch=='T' || ch=='t') return 'A';
    if(ch=='C' || ch=='c') return 'G';
    if(ch=='G' || ch=='g') return 'C';
    return 'N';
}

char *dna_rev(char *sr,const char *s,int len){
    int i;
    for(i=0;i<len;i++){
        sr[i]=dna_rev_char(s[len-i-1]);
    }
    sr[i]=0;
    return sr;
}

int readLine(FILE *file){
    char ch;
    while(fscanf(file,"%c",&ch)==1) if(ch=='\n') return 1;
    return 0;
}

char str[GENE_LEN];
char buffer[9182];

int main(void){
    int count=0;
    MYSQL *my_conn;
    MYSQL_RES *result;
    MYSQL_ROW sql_row;

    my_conn=mysql_init(NULL);
    if(mysql_real_connect(my_conn,"127.0.0.1","root","zy19930108","db",3306,NULL,0)){
        FILE *fout=fopen("d:/out.sql","w");
        char path[]="C:/Users/ZhaoYi/Desktop/NNY-Database/Database/";
        strcpy(buffer,path);
        strcat(buffer,"req.txt");
        FILE *req=fopen(buffer,"r");

        char name[100],pam[20],dp[512];
        char cname[100],fp[1024];
        int Chr_No;
        while(fscanf(req,"%s\t%s\t%s",name,pam,dp)==3){
            int req_pam_len=strlen(pam);
            strcpy(buffer,path);
            strcat(buffer,dp);
            strcat(buffer,"list.txt");
            FILE *list=fopen(buffer,"r");

            while(fscanf(list,"%s\t%s",cname,fp)==2){
                char nt[LEN+1],ch;

                fscanf(list,"%s",buffer);
                fscanf(list,"%s",buffer);
                sprintf(buffer,"select Chr_No from Table_Specie as s join Table_chromosome as c where s.Sno=c.Sno and SName='%s' and Chr_Name='%s';",name,cname);
                int res=mysql_query(my_conn,buffer);
                if(!res){
                    result=mysql_store_result(my_conn);
                    if(result){
                        sql_row=mysql_fetch_row(result);
                        sscanf(sql_row[0],"%d",&Chr_No);
                    }
                }

                strcpy(buffer,path);
                strcat(buffer,dp);
                strcat(buffer,fp);
                FILE *ff=fopen(buffer,"r");
                readLine(ff);
                int i=1,j;
                int len;
                while(fscanf(ff,"%c",&ch)==1){
                    if(ch==10) continue;
                    str[i++]=ch;
                }
                len=i;
                fclose(ff);

                for(i=LEN;i<len-req_pam_len;i++){       // All possible gRNAs, +direction
                    if(check_pam(str+i,pam)){
                        for(j=0;j<LEN;j++) nt[j]=(str+i-LEN)[j];
                        nt[j]=0;
                        sprintf(buffer,"\"%d\",\"%d\",\"%d\",\"%c\",\"%s\",\"%s\"",Chr_No,i,i+LEN-1,'+',nt,pam);
                        fprintf(fout,"%s\n",buffer);
                        count++;
                    }
                }
                char req_pam_rev[PAM_LEN];
                dna_rev(req_pam_rev,pam,req_pam_len);
                for(i=0;i<len-LEN-req_pam_len;i++){     // All possible gRNAs, -direction
                    if(check_pam(str+i,req_pam_rev)){
                        for(j=0;j<LEN;j++) nt[j]=dna_rev_char((str+i+req_pam_len)[LEN-j-1]);
                        nt[j]=0;
                        sprintf(buffer,"\"%d\",\"%d\",\"%d\",\"%c\",\"%s\",\"%s\"",Chr_No,i+req_pam_len-LEN,i+req_pam_len-1,'-',nt,pam);
                        fprintf(fout,"%s\n",buffer);
                        count++;
                    }
                }
            }
        }
        int res2=mysql_query(my_conn,"load data infile 'd:/out.sql' ignore into table Table_sgrna fields terminated by ',' enclosed by '\"' lines terminated by '\r\n' (Chr_No, sgrna_start, sgrna_end, sgrna_strand, sgrna_seq, sgrna_PAM);");
        if(res2!=0) fprintf(stderr,"Load error %d:%s\n",mysql_errno(my_conn),mysql_error(my_conn));
        fclose(fout);
    }
    printf("%d rows will be affected\n",count);
    mysql_free_result(result);
    mysql_close(my_conn);

    return 0;
}
