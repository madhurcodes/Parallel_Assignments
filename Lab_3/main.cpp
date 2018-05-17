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

#include<bitset>

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
    // PAD
    int pad_num_rows = psize - (maxnumfloats % psize);

    int num_real_cols = num_columns;
    maxnumfloats = maxnumfloats + pad_num_rows;
    int pad_num_cols = psize - (num_columns % psize);
    num_columns = num_columns + pad_num_cols;
    int n_excess_rows = maxnumfloats % psize;
    int n_my_rows;
    if (prank < n_excess_rows){
        n_my_rows = (maxnumfloats/psize) + 1;
    }
    else{
        n_my_rows = (maxnumfloats/psize);
    }
    int n_excess_columns = num_columns % psize;
    int n_my_columns;
    if (prank < n_excess_columns){
        n_my_columns = (num_columns/psize) + 1;
    }
    else{
        n_my_columns = (num_columns/psize);
    }




    vector <int> col_to_proc_rank(num_columns);
    vector <int> col_to_proc_off(num_columns);
    // cout <<"test"<<endl;
    int j = 0;
    int k = 0;
    for (i=0;i<num_columns;i++){

        col_to_proc_rank[i] = j;
        col_to_proc_off[i] = k;
        // cout<<"i is "<<i<<" col_proc_r is " <<j<<" off is "<<k<<endl;
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

    // vector<vector<pair<float,char *> > > mycolumns;
    vector<vector<pair<float,pair<char *,int> > > > mycolumns;

    int my_columns_string_lengths[n_my_columns];
    char **val;
    // cout<< "RANK is "<<prank<<endl;
    for(i=0;i<n_my_columns;i++){
        
        char fname[80];
        char num[10];
        // if (prank<n_excess_columns){
        //     sprintf(num,"%d",prank*((num_columns/psize) +1) + i + 1 );
        //     // cout<< "RANK is "<<prank<<" I is "<<i<<" "<<prank*((num_columns/psize) +1) + i + 1 <<endl;
        // }
        // else{
        //     sprintf(num,"%d",prank*((num_columns/psize)+1) - (prank-n_excess_columns)  + i + 1 );
        //     // cout<< "RANK is "<<prank<<" I is "<<i<<" "<<prank*((num_columns/psize)+1) - (prank-n_excess_columns)  + i + 1 <<endl;
        // }
        int inp_i;
        if(prank*(num_columns/psize)  + i + 1 > num_real_cols ){

            int num_floats = 0;
            my_columns_string_lengths[i] = 0;
            vector<pair<float,pair<char *,int> > > column;
            int j = 0;
            float inp_f;
            inp_i = my_columns_string_lengths[i] ;
            // char val[maxnumfloats][maxlengthstring+1] = {0};
            // memset(val,0,sizeof(val[0][0])*maxnumfloats*(maxlengthstring+1));
            val = (char**) calloc(maxnumfloats, sizeof(char*));;
            for (j= 0; j< maxnumfloats; j++ )
            {
            val[j] = (char*) calloc((maxlengthstring+1), sizeof(char));
            }
            for(j=0;j<maxnumfloats;j++){
                inp_f = numeric_limits<float>::max();
                pair<float,pair<char *,int> > kv(inp_f,make_pair(val[j], inp_i ));
                column.push_back(kv);
            }
            // cout<<"alloc as "<<inp_i<<endl;
            free(val);
            mycolumns.push_back(column);
        }
        else{
            sprintf(num,"%d",prank*(num_columns/psize)  + i + 1 );

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

            vector<pair<float,pair<char *,int> >  > column;
            int j = 0;
            float inp_f;
            inp_i = my_columns_string_lengths[i] ;
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
                pair<float,pair<char *,int> >  kv(inp_f,make_pair(val[j], inp_i ));
                // cout<<"alloc as "<<inp_i<<endl;
                column.push_back(kv);
            }
            // free(val);
            fclose(file);
            free(val);
            mycolumns.push_back(column);
        }
    }


    vector<vector<pair<float,pair<char *,int> > > > myrows;

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
    char *csendbfr ;
    char *crecvbfr;
    int sendcounts[psize], sdispls[psize],  recvcounts[psize],rdispls[psize];
    int csendcounts[psize], csdispls[psize],  crecvcounts[psize],crdispls[psize]; 

        // cout<<"========================================"<<endl;
        // cout<<" alloc cols   "<<endl;
        // trank = 0;
        // while (trank < psize) {
        //     if (prank == trank) {
        //         printf ("col printed by rank: %d\n", trank);
        //         for (i=0;i<mycolumns.size();i++){
        //             for(k=0;k<mycolumns[i].size();k++){
        //                 cout<<" |k = "<<mycolumns[i][k].first<<" v = "<<mycolumns[i][k].second.first<<" s is " <<mycolumns[i][k].second.second;
        //             }
        //             cout<<endl<<"------------------------------------------------"<<endl;
        //         }
        //         fflush (stdout);
        //     }
        //     trank ++;
        //     MPI_Barrier (MPI_COMM_WORLD);
        // }
        // cout<<"========================================"<<endl;


    while (!(DONE==1)){

        //columns are in vectors row is empty so we first transfer stuff to rows

        // sendbfr = (float *) calloc(sizeof(float)*maxnumfloats*n_my_columns);
        // recvbfr = (float *) calloc(sizeof(float)*num_columns*n_my_rows);
        // csendbfr = (char *) calloc((maxlengthstring+1)*sizeof(char)*maxnumfloats*n_my_columns) ;
        // crecvbfr = (char *) calloc((maxlengthstring+1)*sizeof(char)*num_columns*n_my_rows);
        MPI_Barrier(MPI_COMM_WORLD);
        sendbfr = (float *) calloc(maxnumfloats*n_my_columns,sizeof(float));
        recvbfr = (float *) calloc(num_columns*n_my_rows,sizeof(float));
        csendbfr = (char *) calloc((maxlengthstring+1 + 4)*maxnumfloats*n_my_columns,sizeof(char)) ;
        crecvbfr = (char *) calloc((maxlengthstring+1 + 4 )*num_columns*n_my_rows,sizeof(char));

        // first putting first elements of all columns in this proc then second elements and so on
        // cout<<"========================================"<<endl;
        // cout<<" alloc ** cols   "<<endl;
        // trank = 0;
        // while (trank < psize) {
        //     if (prank == trank) {
        //         printf ("col printed by rank: %d\n", trank);
        //         for (i=0;i<mycolumns.size();i++){
        //             for(k=0;k<mycolumns[i].size();k++){
        //                 cout<<" |k = "<<mycolumns[i][k].first<<" v = "<<mycolumns[i][k].second.first<<" s is " <<mycolumns[i][k].second.second;
        //             }
        //             cout<<endl<<"------------------------------------------------"<<endl;
        //         }
        //         fflush (stdout);
        //     }
        //     trank ++;
        //     MPI_Barrier (MPI_COMM_WORLD);
        // }
        // cout<<"========================================"<<endl;


        for(i=0;i<mycolumns.size();i++){
            for(j=0;j<mycolumns[i].size();j++){
                float temp_f = mycolumns[i][j].first;
                sendbfr[j*n_my_columns+i] = temp_f;
                for(k=0;k<maxlengthstring+1;k++){
                    char temp_c = mycolumns[i][j].second.first[k];
                    csendbfr[(j*n_my_columns+i) * (maxlengthstring+1+4) + k] = temp_c;
                }
                int temp_i = mycolumns[i][j].second.second;
                char *seri = reinterpret_cast<char*> (&temp_i);
                int lv = 0;
                for(k=maxlengthstring+1;k<maxlengthstring+1+ 4;k++){
                    csendbfr[(j*n_my_columns+i) * (maxlengthstring+1+4) + k] = seri[lv];
                    // cout<<"send as "<<  <<endl;
                    lv++;
                }
                // cout<<"i cud haz "<< mycolumns[i][j].second.second<< endl;
                // cout<<" i send ax "<< *(reinterpret_cast<int *> (&csendbfr[(j*n_my_columns+i) * (maxlengthstring+1) + maxlengthstring+1])) <<endl;

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
            csendcounts[i] = n_my_columns*(row_counts[i]  )*(maxlengthstring+1 + 4);
            crecvcounts[i] = column_counts[i]*n_my_rows*(maxlengthstring+1 + 4);
            // cout<<"SSS is "<<sizeof(int)<<endl;

            if(i==0){
                csdispls[i]  = 0;
                crdispls[i] = 0 ;
            }
            else{
                // csdispls[i] = csdispls[i-1] += sizeof(char)*n_my_columns*temp_num_rows*(maxlengthstring+1);
                // crdispls[i] = crdispls[i-1] += sizeof(char)*n_my_rows*temp_num_columns*(maxlengthstring+1);
                // csdispls[i] = csdispls[i-1] + n_my_columns*row_counts[i] *(maxlengthstring+1);
                // crdispls[i] = crdispls[i-1] + n_my_rows*column_counts[i]*(maxlengthstring+1);
                csdispls[i] = csdispls[i-1] + n_my_columns*row_counts[i-1] *(maxlengthstring+1 + 4);
                crdispls[i] = crdispls[i-1] + n_my_rows*column_counts[i-1]*(maxlengthstring+1+ 4 );


            }

        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Alltoallv(sendbfr, sendcounts, sdispls, MPI_FLOAT, recvbfr, recvcounts, rdispls, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Alltoallv(csendbfr, csendcounts, csdispls, MPI_CHAR, crecvbfr, crecvcounts, crdispls, MPI_CHAR, MPI_COMM_WORLD);
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        for(i=0;i<n_my_rows;i++){
            vector<pair<float,pair<char *,int> > > row;
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
            int inp_i;

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
                    val[j][k] = crecvbfr[ crdispls[proc_rank] + (maxlengthstring+1+4)*(i*column_counts[proc_rank] + col_to_proc_off[j] ) + k ];                    
                    // if(proc_rank != 0){
                    // val[j][k] = crecvbfr[ crdispls[proc_rank-1] + (maxlengthstring+1)*(i*column_counts[proc_rank] + col_to_proc_off[j] ) + k ];
                    // }
                    // else{
                    // val[j][k] = crecvbfr[ (maxlengthstring+1)*(i*column_counts[proc_rank] + col_to_proc_off[j] ) + k ];
                    // }
                }
                inp_i =  *( reinterpret_cast<int *> (&crecvbfr[ crdispls[proc_rank] + (maxlengthstring+1+4)*(i*column_counts[proc_rank] + col_to_proc_off[j] ) + maxlengthstring+1 ] ));
                // cout <<"recv asm "<<inp_i<<endl;
                pair<float,pair<char *,int> > kv(inp_f,make_pair(val[j],inp_i)  );
                row.push_back(kv);
                // cout<<inp_f<<val[j]<<endl;
            }
            free(val);
            myrows.push_back(row);
        }
        free(sendbfr);
        free(csendbfr);
        free(recvbfr);
        free(crecvbfr);


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
                free(mycolumns[i][j].second.first);
            }
            mycolumns[i] = vector<pair<float, pair<char *,int> > > ();
        }
        mycolumns = vector<vector<pair<float, pair<char *,int> > > > ();

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

        // MPI_Barrier (MPI_COMM_WORLD);
        // cout<<"========================================"<<endl;
        // cout<<"ROWS   "<<endl;
        // trank = 0;
        // while (trank < psize) {
        //     if (prank == trank) {
        //         printf ("Row printed by rank: %d\n", trank);
        //         for (i=0;i<myrows.size();i++){
        //             for(k=0;k<myrows[i].size();k++){
        //                 cout<<" |k = "<<myrows[i][k].first<<"  v = "<<myrows[i][k].second.first<<" len = "<<myrows[i][k].second.second;
        //             }
        //             cout<<endl<<"------------------------------------------------"<<endl;
        //         }
        //         fflush (stdout);
        //     }
        //     trank ++;    
        //     MPI_Barrier (MPI_COMM_WORLD);
        // }
        // // cout<<"========================================"<<endl;


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
        // csendbfr = (char *) calloc((maxlengthstring+1)*sizeof(char)*num_columns*n_my_rows) ;
        // crecvbfr = (char *) calloc((maxlengthstring+1)*sizeof(char)*maxnumfloats*n_my_columns);
        
        sendbfr = (float *) calloc(num_columns*n_my_rows,sizeof(float));
        recvbfr = (float *) calloc(maxnumfloats*n_my_columns,sizeof(float));
        csendbfr = (char *) calloc((maxlengthstring+1+4)*num_columns*n_my_rows,sizeof(char)) ;
        crecvbfr = (char *) calloc((maxlengthstring+1+4 )*maxnumfloats*n_my_columns,sizeof(char));



        // first putting first elements of all columns in this proc then second elements and so on

        for(i=0;i<myrows.size();i++){
            for(j=0;j<myrows[i].size();j++){
                float temp_f = myrows[i][j].first;
                sendbfr[j*n_my_rows+i] = temp_f;
                for(k=0;k<maxlengthstring+1;k++){
                    char temp_c = myrows[i][j].second.first[k];
                    csendbfr[(j*n_my_rows+i) * (maxlengthstring+1+4) + k] = temp_c;
                }

                int temp_i = myrows[i][j].second.second;
                char *seri = reinterpret_cast<char*>( &temp_i);
                int lv = 0;
                for(k=maxlengthstring+1;k<maxlengthstring+1+ 4;k++){
                    csendbfr[(j*n_my_rows+i) * (maxlengthstring+1+4)  + k] = seri[lv];
                    lv++;
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

            csendcounts[i] = column_counts[i]*n_my_rows*(maxlengthstring+1 + 4) ;
            crecvcounts[i] = n_my_columns*(row_counts[i] )*(maxlengthstring+1 + 4);
            if(i==0){
                csdispls[i]  = 0;
                crdispls[i] = 0 ;
            }
            else{
                // csdispls[i] = csdispls[i-1] += sizeof(char)*n_my_columns*temp_num_rows*(maxlengthstring+1);
                // crdispls[i] = crdispls[i-1] += sizeof(char)*n_my_rows*temp_num_columns*(maxlengthstring+1);
                // csdispls[i] = crdispls[i-1] + n_my_rows*column_counts[i]*(maxlengthstring+1);
                // crdispls[i] =  csdispls[i-1] + n_my_columns*row_counts[i]*(maxlengthstring+1);
                csdispls[i] = crdispls[i-1] + n_my_rows*column_counts[i-1]*(maxlengthstring+1 + 4);
                crdispls[i] =  csdispls[i-1] + n_my_columns*row_counts[i-1]*(maxlengthstring+1 + 4);


            }

        }
 
        //  MPI_Finalize();
        // return 0;

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Alltoallv(sendbfr, sendcounts, sdispls, MPI_FLOAT, recvbfr, recvcounts, rdispls, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Alltoallv(csendbfr, csendcounts, csdispls, MPI_CHAR, crecvbfr, crecvcounts, crdispls, MPI_CHAR, MPI_COMM_WORLD);
        
        MPI_Barrier(MPI_COMM_WORLD);

        // MPI_Finalize();

        // return 0;




        for(i=0;i<n_my_columns;i++){
            vector<pair<float,pair<char *,int> > > column;
            float inp_f;
            int inp_i;
            val = (char**) calloc(maxnumfloats, sizeof(char*));
            
            for (j= 0; j< maxnumfloats; j++ )
            {
                val[j] = (char*) calloc((maxlengthstring+1), sizeof(char));
            }

            for(j=0;j<maxnumfloats;j++){
                int proc_rank = row_to_proc_rank[j];
                //-----
                inp_f = recvbfr[ rdispls[proc_rank] + i*row_counts[proc_rank] + row_to_proc_off[j]  ];
                // if(proc_rank!=0){
                //     inp_f = recvbfr[ rdispls[proc_rank-1] + i*row_counts[proc_rank] + row_to_proc_off[j]  ];
                // }
                // else{
                //     inp_f = recvbfr[  i*row_counts[proc_rank] + row_to_proc_off[j]  ];
                // }
                for(k=0;k<(maxlengthstring+1);k++){
                    val[j][k] = crecvbfr[ crdispls[proc_rank] + (maxlengthstring+1+4)*(i*row_counts[proc_rank] + row_to_proc_off[j] ) + k ];
                    // if(proc_rank != 0){
                    // val[j][k] = crecvbfr[ crdispls[proc_rank-1] + (maxlengthstring+1)*(i*row_counts[proc_rank] + row_to_proc_off[j] ) + k ];
                    // }
                    // else{
                    // val[j][k] = crecvbfr[ (maxlengthstring+1)*(i*row_counts[proc_rank] + row_to_proc_off[j] ) + k ];
                    // }
                }
                inp_i =  *(reinterpret_cast<int *> ( &crecvbfr[ crdispls[proc_rank] + (maxlengthstring+1+ 4)*(i*row_counts[proc_rank] + row_to_proc_off[j] )  + maxlengthstring+1 ]) );
                pair<float,pair<char *,int> > kv(inp_f,make_pair(val[j],inp_i)  );
                // pair<float, pair<char *,int> > kv(inp_f,val[j]);
                column.push_back(kv);
            }
            free(val);
            mycolumns.push_back(column);
        }
        free(sendbfr);
        free(csendbfr);
        free(recvbfr);
        free(crecvbfr);

        // myrows.resize(0);
        for (i=0;i<myrows.size();i++){
            for(j=0;j<myrows[i].size();j++){
                free(myrows[i][j].second.first);
            }
            myrows[i] = vector<pair<float,pair<char *,int> > > ();
        }

        myrows = vector<vector<pair<float,pair<char *,int> > > > ();
        

        MPI_Barrier(MPI_COMM_WORLD);
        for(i=0;i<mycolumns.size();i++){
            sort(mycolumns[i].begin(),mycolumns[i].end());
        }
        //-----copy end

    }


    // MPI_Barrier (MPI_COMM_WORLD);
    // cout<<"========================================"<<endl;
    // cout<<"ROWS   "<<endl;
    // trank = 0;
    // while (trank < psize) {
    //     if (prank == trank) {
    //         printf ("Row printed by rank: %d\n", trank);
    //         for (i=0;i<myrows.size();i++){
    //             for(k=0;k<myrows[i].size();k++){
    //                 cout<<" |k = "<<myrows[i][k].first<<"  v = "<<myrows[i][k].second.first<<" len = "<<myrows[i][k].second.second;
    //             }
    //             cout<<endl<<"------------------------------------------------"<<endl;
    //         }
    //         fflush (stdout);
    //     }
    //     trank ++;    
    //     MPI_Barrier (MPI_COMM_WORLD);
    // }
    // cout<<"========================================"<<endl;




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
            if (prank == 0){
                ofstream clear_f(fname, ios::out | ios::trunc);
                clear_f.close();
            }
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

                        my_output.write(myrows[i][j].second.first,(sizeof(char))*myrows[i][j].second.second);
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
                free(myrows[i][j].second.first);
            }
            myrows[i] = vector<pair<float,pair<char *,int> > > ();
        }
    myrows = vector<vector<pair<float,pair<char *,int> > > > ();
    // cout<<"========================================"<<endl;

    // cout<<"My Rank is "<<prank<<endl;
    MPI_Finalize();
    return 0;
}