#include <stdio.h>  /* fopen */
#include <stdlib.h>  /* malloc */
#include "include/cjson.h"
#include "include/hhmpcdynmem.h"

cJSON* hhmpc_dynmem_get_data (char *fname)
{
    cJSON *data;
    char *fdata;
    FILE *f=fopen(fname,"rb");
    if (NULL == f) {
        printf("ERROR: could not open file %s \n", fname);
        return NULL;
    }
    
    fseek(f,0,SEEK_END);
    long len=ftell(f);
    fseek(f,0,SEEK_SET);
    fdata=(char*)malloc(len+1);
    if (NULL == fdata) {
        printf("ERROR: could not allocate memory for file %s \n", fname);
        return NULL;
    }
    
    fread(fdata,1,len,f);
    fclose(f);
    data = cJSON_Parse(fdata);
    free(fdata);
    if (NULL == data) {
        printf("ERROR: could not parse file %s \n", fname);
        return NULL;
    }
    
    return data;
}
