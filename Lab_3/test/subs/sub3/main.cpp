#include <mpi.h>
#include <iostream>
#include <vector>
#include <algorithm> 
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include<limits>
#include <cstring>
#include <fstream>

using namespace std;

int main(int argc , char * argv[]){


    MPI_Init(&argc, &argv);
    int DEBUG_ITERATION = 0;
    int num_columns = atoi(argv[1]);
    char* base_filename = argv[2];

    int trank = 0;

    // char fname[80];
    // char num[10];
    // sprintf(num,"%d",num_columns);
    // strcpy(fname,base_filename);
    // strcat(fname,num);
    // cout<<"fname si "<<fname<<endl;
    int prank; MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    int psize; MPI_Comm_size(MPI_COMM_WORLD, &psize);


    int n_excess_columns = num_columns % psize;
    int n_my_columns;
    if (prank < n_excess_columns){
        n_my_columns = (num_columns/psize) + 1;
    }
    else{
        n_my_columns = (num_columns/psize);
    }

    int maxlengthstring = -1; 
    int maxnumfloats = -1;
    int i = 0;
    for(i=0;i<num_columns;i++){
        char fname[80];
        char num[10];
        sprintf(num,"%d", i + 1 );
        strcpy(fname,base_filename);
        strcat(fname,num);
        FILE *file = fopen(fname, "rb");
        fseek(file, 0, SEEK_END);
        int csize = ftell(file);
        fseek(file, 0, SEEK_SET);
        int length_of_string;

        fread((void*)(&length_of_string), sizeof(length_of_string), 1, file);
        
        int num_floats = (csize-4)/ (length_of_string + 4);
        if (length_of_string>maxlengthstring){
            maxlengthstring = length_of_string;
        }
        if(num_floats>maxnumfloats){
            maxnumfloats = num_floats;
        }

        fclose(file);
        
    }
    int n_excess_rows = maxnumfloats % psize;
    int n_my_rows;
    if (prank < n_excess_rows){
        n_my_rows = (maxnumfloats/psize) + 1;
    }
    else{
        n_my_rows = (maxnumfloats/psize);
    }




    vector <int> col_to_proc_rank(num_columns);
    vector <int> col_to_proc_off(num_columns);
    
    int j = 0;
    int k = 0;
    for (i=0;i<num_columns;i++){

        col_to_proc_rank[i] = j;
        col_to_proc_off[i] = k;
        if (j<n_excess_columns){
            if(k>=(num_columns/psize)){
                j+=1;
                k = 0;
            }
            else{
                k++;
            }
        }
        else{
            if(k>=(num_columns/psize)-1){
                j+=1;
                k = 0;
            }
            else{
                k++;
            }

        }
    }


    vector <int> row_to_proc_rank(maxnumfloats);
    vector <int> row_to_proc_off(maxnumfloats);
    j = 0;
    k = 0;
    for (i=0;i<maxnumfloats;i++){

        row_to_proc_rank[i] = j;
        row_to_proc_off[i] = k;
        if (j<n_excess_rows){
            if(k>=(maxnumfloats/psize)){
                j+=1;
                k = 0;
            }
            else{
                k++;
            }
        }
        else{
            if(k>=(maxnumfloats/psize)-1){
                j+=1;
                k = 0;
            }
            else{
                k++;
            }

        }
    }



    vector<vector<pair<float,char *> > > mycolumns;

    int my_columns_string_lengths[n_my_columns];
    char ** val;
    for(i=0;i<n_my_columns;i++){
        char fname[80];
        char num[10];
        if (prank<n_excess_columns){
            sprintf(num,"%d",prank*((num_columns/psize) +1) + i + 1 );
        }
        else{
            sprintf(num,"%d",prank*((num_columns/psize)+1) - (prank-n_excess_columns)  + i + 1 );
        }
        
        strcpy(fname,base_filename);
        strcat(fname,num);
        FILE *file = fopen(fname, "rb");
        fseek(file, 0, SEEK_END);
        int csize = ftell(file);
        fseek(file, 0, SEEK_SET);
        int length_of_string;

        fread((void*)(&length_of_string), sizeof(length_of_string), 1, file);
        
        int num_floats = (csize-4)/ (length_of_string + 4);
        my_columns_string_lengths[i] = length_of_string;

        vector<pair<float,char *> > column;
        int j = 0;
        float inp_f;
        // char val[maxnumfloats][maxlengthstring+1] = {0};
        // memset(val,0,sizeof(val[0][0])*maxnumfloats*(maxlengthstring+1));
        val = (char**) calloc(maxnumfloats, sizeof(char*));;
        for (j= 0; j< maxnumfloats; j++ )
        {
         val[j] = (char*) calloc((maxlengthstring+1), sizeof(char));
        }
        for(j=0;j<maxnumfloats;j++){
            if (j<num_floats){
                fread((void*)(&inp_f), sizeof(inp_f), 1, file);
                fread((void*)( val[j] ), length_of_string, 1, file);
            }
            else{
                inp_f = numeric_limits<float>::max();
            }
            // cout<<val[i]<<endl;
            pair<float,char *> kv(inp_f,val[j]);
            column.push_back(kv);
        }
        // free(val);
        fclose(file);
        free(val);
        mycolumns.push_back(column);
    }




    // for (i=0;i<mycolumns.size();i++){
    //     string showMsgStr(mycolumns[i][0].second, mycolumns[i][0].second + my_columns_string_lengths[0]);
    //     cout<<mycolumns[0][0].first<<"  " << showMsgStr <<endl;
    //     cout<<"++++++++++++++++++++++++++++++++++"<<endl;
    // }

    vector<vector<pair<float,char *> > > myrows;

    int column_counts[psize];
    int row_counts[psize];

    int string_length_disps[psize];
    for(i=0;i<psize;i++){
        int temp_num_rows; // num rows proc i will have
        int temp_num_columns; // num cols proc i will have
        if(i<n_excess_rows){
            temp_num_rows = (maxnumfloats/psize) + 1;
        }
        else{
            temp_num_rows = (maxnumfloats/psize);
        }
        if(i<n_excess_columns){
            temp_num_columns = (num_columns/psize) + 1;
        }
        else{
            temp_num_columns = (num_columns/psize);
        }
        if(i==0){
            string_length_disps[i] = 0;
        }
        else{
            string_length_disps[i] = string_length_disps[i-1] + column_counts[i-1];
        }
        column_counts[i] = temp_num_columns;
        row_counts[i] = temp_num_rows;
    }

    int DONE = 0;

    int global_string_length[num_columns];

    MPI_Allgatherv(my_columns_string_lengths,n_my_columns,MPI_INT,global_string_length,column_counts,string_length_disps,MPI_INT,MPI_COMM_WORLD);
    
    // cout<< "SL's are "<<endl;
    // for(i=0;i<num_columns;i++){
    //     cout<<" "<<global_string_length[i]<<" ";
    // }
    // cout<<"++++++++"<<endl;


    float *sendbfr ;
    float *recvbfr ;
    char *sendcbfr ;
    char *recvcbfr;
    int sendcounts[psize], sdispls[psize],  recvcounts[psize],rdispls[psize];
    int csendcounts[psize], csdispls[psize],  crecvcounts[psize],crdispls[psize]; 


    while (!(DONE==1)){

        //columns are in vectors row is empty so we first transfer stuff to rows

        // sendbfr = (float *) calloc(sizeof(float)*maxnumfloats*n_my_columns);
        // recvbfr = (float *) calloc(sizeof(float)*num_columns*n_my_rows);
        // sendcbfr = (char *) calloc((maxlengthstring+1)*sizeof(char)*maxnumfloats*n_my_columns) ;
        // recvcbfr = (char *) calloc((maxlengthstring+1)*sizeof(char)*num_columns*n_my_rows);
        MPI_Barrier(MPI_COMM_WORLD);
        sendbfr = (float *) calloc(maxnumfloats*n_my_columns,sizeof(float));
        recvbfr = (float *) calloc(num_columns*n_my_rows,sizeof(float));
        sendcbfr = (char *) calloc((maxlengthstring+1)*maxnumfloats*n_my_columns,sizeof(char)) ;
        recvcbfr = (char *) calloc((maxlengthstring+1)*num_columns*n_my_rows,sizeof(char));

        // first putting first elements of all columns in this proc then second elements and so on

        for(i=0;i<mycolumns.size();i++){
            for(j=0;j<mycolumns[i].size();j++){
                float temp_f = mycolumns[i][j].first;
                sendbfr[j*n_my_columns+i] = temp_f;
                for(k=0;k<maxlengthstring+1;k++){
                    char temp_c = mycolumns[i][j].second[k];
                    sendcbfr[(j*n_my_columns+i) * (maxlengthstring+1) + k] = temp_c;
                }
            }
        }
        // int sendcounts[psize], sdispls[psize],  recvcounts[psize],rdispls[psize];
        // int csendcounts[psize], csdispls[psize],  crecvcounts[psize],crdispls[psize]; 



        // how many cols/roes proc i will have

        // cout<<"row counts is "<<row_counts[0]<< " "<<row_counts[1]<<endl;
        // cout<<"col counts is "<<column_counts[0]<< " "<<column_counts[1]<<endl;

        for(i=0;i<psize;i++){

            sendcounts[i] = n_my_columns*(row_counts[i]  );
            recvcounts[i] = column_counts[i]*n_my_rows;
            if(i==0){
                sdispls[i]  = 0;
                rdispls[i] = 0 ;
            }
            else{
                // sdispls[i] = sdispls[i-1] += sizeof(float)*n_my_columns*temp_num_rows;
                // rdispls[i] = rdispls[i-1] += sizeof(float)*n_my_rows*temp_num_columns;
                // sdispls[i] = sdispls[i-1] + n_my_columns*row_counts[i] ;
                // rdispls[i] = rdispls[i-1] + n_my_rows*column_counts[i];
                sdispls[i] = sdispls[i-1] + n_my_columns*row_counts[i-1] ;
                rdispls[i] = rdispls[i-1] + n_my_rows*column_counts[i-1];

            }

            // csendcounts[i] = n_my_columns*(row_counts[i]  )*(maxlengthstring+1);
            // crecvcounts[i] = column_counts[i]*n_my_rows*(maxlengthstring+1);
            csendcounts[i] = n_my_columns*(row_counts[i]  )*(maxlengthstring+1);
            crecvcounts[i] = column_counts[i]*n_my_rows*(maxlengthstring+1);

            if(i==0){
                csdispls[i]  = 0;
                crdispls[i] = 0 ;
            }
            else{
                // csdispls[i] = csdispls[i-1] += sizeof(char)*n_my_columns*temp_num_rows*(maxlengthstring+1);
                // crdispls[i] = crdispls[i-1] += sizeof(char)*n_my_rows*temp_num_columns*(maxlengthstring+1);
                // csdispls[i] = csdispls[i-1] + n_my_columns*row_counts[i] *(maxlengthstring+1);
                // crdispls[i] = crdispls[i-1] + n_my_rows*column_counts[i]*(maxlengthstring+1);
                csdispls[i] = csdispls[i-1] + n_my_columns*row_counts[i-1] *(maxlengthstring+1);
                crdispls[i] = crdispls[i-1] + n_my_rows*column_counts[i-1]*(maxlengthstring+1);


            }

        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Alltoallv(sendbfr, sendcounts, sdispls, MPI_FLOAT, recvbfr, recvcounts, rdispls, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Alltoallv(sendcbfr, csendcounts, csdispls, MPI_CHAR, recvcbfr, crecvcounts, crdispls, MPI_CHAR, MPI_COMM_WORLD);
        
        MPI_Barrier(MPI_COMM_WORLD);
        for(i=0;i<n_my_rows;i++){
            vector<pair<float,char *> > row;
            int j = 0;
            float inp_f;
            // char val[maxnumfloats][maxlengthstring+1] = {0};
            // memset(val,0,sizeof(val[0][0])*maxnumfloats*(maxlengthstring+1));
            val = (char**) calloc(num_columns, sizeof(char*));
            for (j= 0; j< num_columns; j++ )
            {
            val[j] = (char*) calloc((maxlengthstring+1), sizeof(char));
            }
            // char ** val = (char**) calloc(maxnumfloats*(maxlengthstring+1), sizeof(char));        

            for(j=0;j<num_columns;j++){
                int proc_rank = col_to_proc_rank[j];
                // if(proc_rank!=0){
                //     inp_f = recvbfr[ rdispls[proc_rank-1] + i*column_counts[proc_rank] + col_to_proc_off[j]  ];
                // }
                // else{
                //     inp_f = recvbfr[  i*column_counts[proc_rank] + col_to_proc_off[j]  ];
                // }
                inp_f = recvbfr[ rdispls[proc_rank] + i*column_counts[proc_rank] + col_to_proc_off[j]  ];

                for(k=0;k<(maxlengthstring+1);k++){
                    val[j][k] = recvcbfr[ crdispls[proc_rank] + (maxlengthstring+1)*(i*column_counts[proc_rank] + col_to_proc_off[j] ) + k ];                    
                    // if(proc_rank != 0){
                    // val[j][k] = recvcbfr[ crdispls[proc_rank-1] + (maxlengthstring+1)*(i*column_counts[proc_rank] + col_to_proc_off[j] ) + k ];
                    // }
                    // else{
                    // val[j][k] = recvcbfr[ (maxlengthstring+1)*(i*column_counts[proc_rank] + col_to_proc_off[j] ) + k ];
                    // }
                }
                pair<float,char *> kv(inp_f,val[j]);
                row.push_back(kv);
                // cout<<inp_f<<val[j]<<endl;
            }
            free(val);
            myrows.push_back(row);
        }
        free(sendbfr);
        free(sendcbfr);
        free(recvbfr);
        free(recvcbfr);


        // cout<<"========================================"<<endl;
        // cout<<"cols   "<<endl;
        // trank = 0;
        // while (trank < psize) {
        //     if (prank == trank) {
        //         printf ("col printed by rank: %d\n", trank);
        //         for (i=0;i<mycolumns.size();i++){
        //             for(k=0;k<mycolumns[i].size();k++){
        //                 cout<<" |k = "<<mycolumns[i][k].first<<" v = "<<mycolumns[i][k].second;
        //             }
        //             cout<<endl<<"------------------------------------------------"<<endl;
        //         }
        //         fflush (stdout);
        //     }
        //     trank ++;
        //     MPI_Barrier (MPI_COMM_WORLD);
        // }
        // cout<<"========================================"<<endl;



        // clear columns
        for (i=0;i<mycolumns.size();i++){
            for(j=0;j<mycolumns[i].size();j++){
                free(mycolumns[i][j].second);
            }
            mycolumns[i] = vector<pair<float,char *> > ();
        }
        mycolumns = vector<vector<pair<float,char *> > > ();

        MPI_Barrier(MPI_COMM_WORLD);

        // MPI_Barrier (MPI_COMM_WORLD);
        // cout<<"========================================"<<endl;
        // cout<<"ROWS   "<<endl;
        // trank = 0;
        // while (trank < psize) {
        //     if (prank == trank) {
        //         printf ("Row printed by rank: %d\n", trank);
        //         for (i=0;i<myrows.size();i++){
        //             for(k=0;k<myrows[i].size();k++){
        //                 cout<<" |k = "<<myrows[i][k].first<<" v = "<<myrows[i][k].second;
        //             }
        //             cout<<endl<<"------------------------------------------------"<<endl;
        //         }
        //         fflush (stdout);
        //     }
        //     trank ++;
        //     MPI_Barrier (MPI_COMM_WORLD);
        // }
        // cout<<"========================================"<<endl;


        // CHECK IF SORTED

        
        int  local_sorted = 1;
        for(i=0;i<myrows.size();i++){
            for (j=0;j<myrows[i].size()-1;j++){
                if( ! (myrows[i][j].first<= myrows[i][j+1].first)){
                // cout <<endl<<"My rank = "<<prank<<" and bad array is ";
                // for(k=0;k<myrows[i].size();k++){
                //     cout <<myrows[i][k].first<<" ";
                // }
                // cout <<" and j is "<<j<<endl;
                local_sorted = 0;
                }
            }
            // if(!is_sorted(myrows[i].begin(),myrows[i].end())){
            //     cout <<endl<<"My rank = "<<prank<<" and bad array is ";
            //     for(j=0;j<myrows[i].size();j++){
            //         cout <<myrows[i][j].first<<" ";
            //     }
            //     cout <<endl;
            //     local_sorted = 0;
            // }
        }
        int global_sorted;
        MPI_Allreduce(&local_sorted, &global_sorted, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);


        // cout<<"  iteration "<<DEBUG_ITERATION<<endl;
        DEBUG_ITERATION++;
        // cout<<"global sorted is "<<global_sorted<<endl;
        if (global_sorted){
            DONE = 1; 
            continue;
        }

        // Sort by rows 
        for(i=0;i<myrows.size();i++){
            sort(myrows[i].begin(),myrows[i].end());
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // now transfer to columns




        //-------------

        //rows are in vectors and  mycolumns is empty so we first transfer stuff to cols



        // sendbfr = (float *) calloc(sizeof(float)*num_columns*n_my_rows);
        // recvbfr = (float *) calloc(sizeof(float)*maxnumfloats*n_my_columns);
        // sendcbfr = (char *) calloc((maxlengthstring+1)*sizeof(char)*num_columns*n_my_rows) ;
        // recvcbfr = (char *) calloc((maxlengthstring+1)*sizeof(char)*maxnumfloats*n_my_columns);
        
        sendbfr = (float *) calloc(num_columns*n_my_rows,sizeof(float));
        recvbfr = (float *) calloc(maxnumfloats*n_my_columns,sizeof(float));
        sendcbfr = (char *) calloc((maxlengthstring+1)*num_columns*n_my_rows,sizeof(char)) ;
        recvcbfr = (char *) calloc((maxlengthstring+1)*maxnumfloats*n_my_columns,sizeof(char));



        // first putting first elements of all columns in this proc then second elements and so on

        for(i=0;i<myrows.size();i++){
            for(j=0;j<myrows[i].size();j++){
                float temp_f = myrows[i][j].first;
                sendbfr[j*n_my_rows+i] = temp_f;
                for(k=0;k<maxlengthstring+1;k++){
                    char temp_c = myrows[i][j].second[k];
                    sendcbfr[(j*n_my_rows+i) * (maxlengthstring+1) + k] = temp_c;
                }
            }
        }


        // free(sendcounts);free(sdispls);free(recvcounts);free(rdispls);
        // free(csendcounts);free(csdispls);free(crecvcounts);free(crdispls);
        // sendcounts[psize] = calloc(sizeof(int)*psize);
        // sdispls[psize] = calloc(sizeof(int)*psize);
        // recvcounts[psize] = calloc(sizeof(int)*psize);
        // rdispls[psize] = calloc(sizeof(int)*psize);
        // csendcounts[psize] = calloc(sizeof(int)*psize);
        // csdispls[psize] = calloc(sizeof(int)*psize);
        // crecvcounts[psize] = calloc(sizeof(int)*psize);
        // crdispls[psize] = calloc(sizeof(int)*psize); 



        for(i=0;i<psize;i++){
            sendcounts[i] = column_counts[i]*n_my_rows  ;
            recvcounts[i] =  n_my_columns*(row_counts[i]  );
            if(i==0){
                sdispls[i]  = 0;
                rdispls[i] = 0 ;
            }
            else{
                // sdispls[i] = sdispls[i-1] += sizeof(float)*n_my_columns*temp_num_rows;
                // rdispls[i] = rdispls[i-1] += sizeof(float)*n_my_rows*temp_num_columns;
                // sdispls[i] = sdispls[i-1] + n_my_rows*column_counts[i]  ;
                // rdispls[i] = rdispls[i-1] + n_my_columns*row_counts[i];
                sdispls[i] = sdispls[i-1] + n_my_rows*column_counts[i-1]  ;
                rdispls[i] = rdispls[i-1] + n_my_columns*row_counts[i-1];

            }

            csendcounts[i] = column_counts[i]*n_my_rows*(maxlengthstring+1) ;
            crecvcounts[i] = n_my_columns*(row_counts[i] )*(maxlengthstring+1);
            if(i==0){
                csdispls[i]  = 0;
                crdispls[i] = 0 ;
            }
            else{
                // csdispls[i] = csdispls[i-1] += sizeof(char)*n_my_columns*temp_num_rows*(maxlengthstring+1);
                // crdispls[i] = crdispls[i-1] += sizeof(char)*n_my_rows*temp_num_columns*(maxlengthstring+1);
                // csdispls[i] = crdispls[i-1] + n_my_rows*column_counts[i]*(maxlengthstring+1);
                // crdispls[i] =  csdispls[i-1] + n_my_columns*row_counts[i]*(maxlengthstring+1);
                csdispls[i] = crdispls[i-1] + n_my_rows*column_counts[i-1]*(maxlengthstring+1);
                crdispls[i] =  csdispls[i-1] + n_my_columns*row_counts[i-1]*(maxlengthstring+1);


            }

        }
        
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Alltoallv(sendbfr, sendcounts, sdispls, MPI_FLOAT, recvbfr, recvcounts, rdispls, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Alltoallv(sendcbfr, csendcounts, csdispls, MPI_CHAR, recvcbfr, crecvcounts, crdispls, MPI_CHAR, MPI_COMM_WORLD);
        
        MPI_Barrier(MPI_COMM_WORLD);

        for(i=0;i<n_my_columns;i++){
            vector<pair<float,char *> > column;
            float inp_f;
//mark  
            val = (char**) calloc(maxnumfloats, sizeof(char*));
            
            for (j= 0; j< maxnumfloats; j++ )
            {
            val[j] = (char*) calloc((maxlengthstring+1), sizeof(char));
            }

            for(j=0;j<maxnumfloats;j++){
                int proc_rank = row_to_proc_rank[j];
                inp_f = recvbfr[ rdispls[proc_rank] + i*row_counts[proc_rank] + row_to_proc_off[j]  ];
                // if(proc_rank!=0){
                //     inp_f = recvbfr[ rdispls[proc_rank-1] + i*row_counts[proc_rank] + row_to_proc_off[j]  ];
                // }
                // else{
                //     inp_f = recvbfr[  i*row_counts[proc_rank] + row_to_proc_off[j]  ];
                // }
                for(k=0;k<(maxlengthstring+1);k++){
                    val[j][k] = recvcbfr[ crdispls[proc_rank] + (maxlengthstring+1)*(i*row_counts[proc_rank] + row_to_proc_off[j] ) + k ];
                    // if(proc_rank != 0){
                    // val[j][k] = recvcbfr[ crdispls[proc_rank-1] + (maxlengthstring+1)*(i*row_counts[proc_rank] + row_to_proc_off[j] ) + k ];
                    // }
                    // else{
                    // val[j][k] = recvcbfr[ (maxlengthstring+1)*(i*row_counts[proc_rank] + row_to_proc_off[j] ) + k ];
                    // }
                }
                pair<float,char *> kv(inp_f,val[j]);
                column.push_back(kv);
            }
            free(val);
            mycolumns.push_back(column);
        }
        free(sendbfr);
        free(sendcbfr);
        free(recvbfr);
        free(recvcbfr);

        // myrows.resize(0);
        for (i=0;i<myrows.size();i++){
            for(j=0;j<myrows[i].size();j++){
                free(myrows[i][j].second);
            }
            myrows[i] = vector<pair<float,char *> > ();
        }

        myrows = vector<vector<pair<float,char *> > > ();
        

        MPI_Barrier(MPI_COMM_WORLD);
        for(i=0;i<mycolumns.size();i++){
            sort(mycolumns[i].begin(),mycolumns[i].end());
        }
        //-----copy end
    }


    MPI_Barrier (MPI_COMM_WORLD);
    // cout<<"========================================"<<endl;
    // cout<<"ROWS   "<<endl;
    // trank = 0;
    // while (trank < psize) {
    //     if (prank == trank) {
    //         printf ("Row printed by rank: %d\n", trank);
    //         for (i=0;i<myrows.size();i++){
    //             for(k=0;k<myrows[i].size();k++){
    //                 cout<<" |k = "<<myrows[i][k].first<<" v = "<<myrows[i][k].second;
    //             }
    //             cout<<endl<<"------------------------------------------------"<<endl;
    //         }
    //         fflush (stdout);
    //     }
    //     trank ++;    
    //     MPI_Barrier (MPI_COMM_WORLD);
    // }
    // // cout<<"========================================"<<endl;




    // after exiting myrows is full so we just write it 
    MPI_Barrier (MPI_COMM_WORLD);
    trank = 0;



    char fname[80];
    char num[15];
    sprintf(num,"%d", 0);
    strcpy(fname,base_filename);
    strcat(fname,num);


    while (trank < psize) {
        if (prank == trank) {
            ofstream my_output(fname,  ios::out | ios::binary | ios::app); 
            // my_output.write((char *)(&trank ), sizeof(trank));
            for (i=0;i<myrows.size();i++){
                for(j=0;j<myrows[i].size();j++){
                    if (myrows[i][j].first == numeric_limits<float>::max()){
                        // cout << "i skipped!!!!"<<endl;
                        continue;
                    }
                    else{
                        my_output.write((char *)(&(myrows[i][j].first) ), sizeof(myrows[i][j].first));

                        my_output.write(myrows[i][j].second,(sizeof(char))*global_string_length[j]);
                        // my_output.write(myrows[i][j].second,5);
                    }


                    // cout <<my_columns_string_lengths[j]<<endl;
                    // cout<<" |k = "<<myrows[i][j].first<<" v = "<<myrows[i][j].second;
                }
            }
            my_output.close();
        }
        trank ++;    
        MPI_Barrier (MPI_COMM_WORLD);
    }

    for (i=0;i<myrows.size();i++){
        for(j=0;j<myrows[i].size();j++){
                free(myrows[i][j].second);
            }
            myrows[i] = vector<pair<float,char *> > ();
        }
    myrows = vector<vector<pair<float,char *> > > ();
    // cout<<"========================================"<<endl;

    // cout<<"My Rank is "<<prank<<endl;
    MPI_Finalize();
    return 0;
}