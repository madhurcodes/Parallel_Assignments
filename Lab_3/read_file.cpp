#include<iostream>
#include<vector>
#include <cstdio>
#include <algorithm>
#include <bitset>
#include<string>

using namespace std;

int main(){

    FILE *file = fopen("col4", "rb");
    fseek(file, 0, SEEK_END);
    int csize = ftell(file);
    fseek(file, 0, SEEK_SET);
    int length_of_string;
    fread((void*)(&length_of_string), sizeof(length_of_string), 1, file);
    int num_floats = (csize-4)/ (length_of_string + 4);

    vector<pair<float,char *> > column;
    int i = 0;
    float num;
    char val[num_floats][length_of_string];
    for(i=0;i<num_floats;i++){
        fread((void*)(&num), sizeof(num), 1, file);
        fread((void*)(&val[i]), length_of_string, 1, file);
        cout<<val[i]<<endl;
        pair<float,char *> kv(num,val[i]);
        column.push_back(kv);
    }
    fclose(file);
    
    cout<<"Size is "<<csize<<" and first is "<<column[0].first<<" num floats "<<column.size()<<endl;   
    cout << "string  is "<<column[0].second <<endl; 

    
    sort(column.begin(), column.end());

    for (i=0; i<column.size(); i++)
    {   
        string showMsgStr(column[i].second, column[i].second + length_of_string);
        cout << column[i].first << " "
             << showMsgStr << endl;
    }
    
    //     vector<int> my_columns_string_lengths;
    // for(i=0;i<n_my_columns;i++){
    //     char fname[80];
    //     char num[10];
    //     if (prank<n_excess_columns){
    //         sprintf(num,"%d",prank*((num_columns/psize) +1) + i + 1 );
    //     }
    //     else{
    //         sprintf(num,"%d",prank*((num_columns/psize)+1) - (prank-n_excess_columns)  + i + 1 );
    //     }
        
    //     strcpy(fname,base_filename);
    //     strcat(fname,num);
    //     FILE *file = fopen(fname, "rb");
    //     fseek(file, 0, SEEK_END);
    //     int csize = ftell(file);
    //     fseek(file, 0, SEEK_SET);
    //     int length_of_string;


    //     fread((void*)(&length_of_string), sizeof(length_of_string), 1, file);
    //     int num_floats = (csize-4)/ (length_of_string + 4);
    //     my_columns_string_lengths.push_back(length_of_string);

    //     vector<pair<float,char *> > column;
    //     int i = 0;
    //     float inp_f;
    //     char val[num_floats][length_of_string];
    //     for(i=0;i<num_floats;i++){
    //         fread((void*)(&inp_f), sizeof(inp_f), 1, file);
    //         fread((void*)(&val[i]), length_of_string, 1, file);
    //         // cout<<val[i]<<endl;
    //         pair<float,char *> kv(inp_f,val[i]);
    //         column.push_back(kv);
    //     }
    //     fclose(file);
    //     mycolumns.push_back(column);

    // }

    return 0;
}