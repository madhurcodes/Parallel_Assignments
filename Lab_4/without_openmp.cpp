#include <mpi.h>

#include <iostream>
#include <vector>
#include <algorithm> 
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <limits>
#include <cstring>
#include <fstream>

#include<bitset>

using namespace std;

int main(int argc , char * argv[]){


    MPI_Init(&argc, &argv);

    char* inp_filename = argv[1];    
    char* out_filename = argv[2];
    int max_rows = atoi(argv[3]);
    int max_columns = atoi(argv[4]);
    
    int prank; MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    int psize; MPI_Comm_size(MPI_COMM_WORLD, &psize);

    FILE *file = fopen(inp_filename, "rb");
    fseek(file, 0, SEEK_END);
    int num_tups = ftell(file)/16;
    fseek(file, 0, SEEK_SET);
    int length_of_string;

    while(max_rows%psize!=0){
        max_rows += 1;
    }

    while(max_columns%psize!=0){
        max_columns += 1;
    }
    // cout<<max_rows<<" "<<max_columns<<endl;
    int num_rows_per_proc = max_rows/psize;
    int num_cols_per_proc = max_columns/psize;
    

    int i,j,k;
    vector<vector< int > > col_indices_in_rows (num_rows_per_proc);
    vector<vector< int > > row_indices_in_cols (num_cols_per_proc);
    

    vector<vector< vector<int> > > list_of_rows (num_rows_per_proc);
    vector<vector< vector<int> > > list_of_cols;
    
    int *sendcounts,*recvcounts,*rdispls,*sdispls;
    sendcounts = (int *) calloc(psize,sizeof(int));
    recvcounts = (int *) calloc(psize,sizeof(int)); 
    rdispls = (int *) calloc(psize,sizeof(int));
    sdispls= (int *) calloc(psize,sizeof(int));

    int ro,co,in;
    float fl;
    for(i=0;i<num_tups;i++){
        fread((void*)(&ro), sizeof(ro), 1, file);
        fread((void*)(&co), sizeof(co), 1, file);
        fread((void*)(&in), sizeof(in), 1, file);
        fread((void*)(&fl), sizeof(fl), 1, file);
        // cout<<"Read "<<ro<<" "<<co<<" "<<fl<<endl;
        // cout<<"Try "<<ro<<" "<<co<<" "<<*reinterpret_cast<float *>(reinterpret_cast<int *>( &fl ))<<endl;
        
        // *reinterpret_cast<float *>( &val );
        // *reinterpret_cast<int *>( &fl );
        if(((co-prank*num_cols_per_proc)>=0) && ((co-prank*num_cols_per_proc)<num_cols_per_proc)){
            recvcounts[ro/num_rows_per_proc] += 4; // div instead of mod
            row_indices_in_cols[co-prank*num_cols_per_proc].push_back(ro);     
        }

        if(((ro-prank*num_rows_per_proc)>=0) && ((ro-prank*num_rows_per_proc)<num_rows_per_proc)){
            col_indices_in_rows[ro-prank*num_rows_per_proc].push_back(co);
            vector <int> var;
            var.push_back(ro);
            var.push_back(co);
            var.push_back(in);
            var.push_back(*reinterpret_cast<int *>( &fl ));

            // cout<<"Var has c= "<<var[1]<<" and pb is ="<<var[3]<<" but fl="<<fl<<endl;
            
            sendcounts[co/num_cols_per_proc] += 4; // 4 ints per entry in array

            list_of_rows[ro-prank*num_rows_per_proc].push_back(var);
        }
    }
    fclose(file);

    for(i=0;i<num_rows_per_proc;i++){
        sort(col_indices_in_rows[i].begin(),col_indices_in_rows[i].end());
    }
    for(i=0;i<num_cols_per_proc;i++){
        sort(row_indices_in_cols[i].begin(),row_indices_in_cols[i].end());
    }

    int sum_send_counts = 0, sum_recv_counts = 0;
    for(i=0;i<psize;i++){
        sdispls[i] = sum_send_counts;
        rdispls[i] = sum_recv_counts;
        sum_send_counts += sendcounts[i];
        sum_recv_counts += recvcounts[i];
    }

    // cout<<"Printing Rows"<<endl;
    // for (i=0;i<list_of_rows.size();i++){
    //     for(k=0;k<list_of_rows[i].size();k++){
    //         cout<<"row - "<<list_of_rows[i][k][0]<<" col - "<<list_of_rows[i][k][1]<< " i - "<<list_of_rows[i][k][2]<< " f - "<<*reinterpret_cast<float *>(&list_of_rows[i][k][3])<<endl;
    //     }
    //     cout<<endl<<"------------------------------------------------"<<endl;
    // }

    int *sendbfr;
    int *recvbfr;
    // int *sendcounts, *sdispls,  *recvcounts,*rdispls;


    int num_it = 0;
    int added_to_s_bfr  = 0;

    for(num_it=0;num_it<4;num_it++){

        //sort rows
        for (i=0;i<list_of_rows.size();i++){
            // argsort one row according to float
            vector <pair<float, pair<int,int> > > tosort;
            for(j=0;j<list_of_rows[i].size();j++){
                tosort.push_back(make_pair(*reinterpret_cast<float *>( &list_of_rows[i][j][3] ), make_pair(list_of_rows[i][j][1],j)));
            }
            sort(tosort.begin(),tosort.end());
            vector <vector<int> > newrow;
            for(j=0;j<list_of_rows[i].size();j++){
                vector<int> newel;
                newel.push_back(list_of_rows[i][tosort[j].second.second][0]);
                newel.push_back(col_indices_in_rows[i][j]);
                newel.push_back(list_of_rows[i][tosort[j].second.second][2]);
                newel.push_back(list_of_rows[i][tosort[j].second.second][3]);
                newrow.push_back(newel);
            }
            for(j=0;j<list_of_rows[i].size();j++){
                for(k=0;k<4;k++){
                    list_of_rows[i][j][k] = newrow[j][k];
                }
            }
            newrow = vector <vector<int> > ();
            tosort = vector <pair<float, pair<int,int> > > ();
        }

        // MPI_Barrier(MPI_COMM_WORLD);
        // int num_els_to_send = 0;
        // for(j=0;j<list_of_rows.size();j++){
        //     num_els_to_send += list_of_rows[j].size();
        // }
        // int sbfr_size = 4*num_els_to_send;


        sendbfr = (int *) calloc(sum_send_counts,sizeof(int));
        recvbfr = (int *) calloc(sum_recv_counts,sizeof(int));
        
        

        added_to_s_bfr  = 0;
        //denotes how far I have travelled along different rows
        vector<int> row_trav(num_rows_per_proc,0);
        for(i=0;i<psize;i++){
            for(j=0;j<list_of_rows.size();j++){
                for(k=row_trav[j];k<list_of_rows[j].size();k++){
                    if((list_of_rows[j][k][1]>=(i*num_cols_per_proc)) && (list_of_rows[j][k][1]<(i+1)*num_cols_per_proc)){
                        sendbfr[added_to_s_bfr] = list_of_rows[j][k][0];
                        sendbfr[added_to_s_bfr+1] = list_of_rows[j][k][1];
                        sendbfr[added_to_s_bfr+2] = list_of_rows[j][k][2];
                        sendbfr[added_to_s_bfr+3] = list_of_rows[j][k][3];
                        added_to_s_bfr += 4;
                    }
                    else{
                        row_trav[j] = k;
                        break;
                    }
                }
            }
        }
        row_trav = vector<int>();

        // for(i=0;i<psize;i++){
        //     for(j=0;j<list_of_rows.size();j++){
        //         for(k=0;k<list_of_rows[j].size();k++){
        //             if((list_of_rows[j][k][1]>=(i*num_cols_per_proc)) && (list_of_rows[j][k][1]<(i+1)*num_cols_per_proc)){
        //                 sendbfr[added_to_s_bfr] = list_of_rows[j][k][0];
        //                 sendbfr[added_to_s_bfr+1] = list_of_rows[j][k][1];
        //                 sendbfr[added_to_s_bfr+2] = list_of_rows[j][k][2];
        //                 sendbfr[added_to_s_bfr+3] = list_of_rows[j][k][3];
        //                 added_to_s_bfr += 4;
        //             }
        //         }
        //     }
        // }
        
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Alltoallv(sendbfr, sendcounts, sdispls, MPI_INT, recvbfr, recvcounts, rdispls, MPI_INT, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        
        
        list_of_rows = vector<vector< vector<int> > > ();
        list_of_cols = vector<vector< vector<int> > > ();
        list_of_cols.resize(num_cols_per_proc);

        // buffer to cols
        for(i=0;i<sum_recv_counts/4;i++){

            int roo = recvbfr[i*4];
            int coo = recvbfr[i*4+1];
            int ino = recvbfr[i*4+2];
            int flo = recvbfr[i*4+3];
            // cout<<"prank = "<< prank<<" i="<<ino <<" roo = "<<roo<<" coo - "<<coo<<endl;
            if(((coo-prank*num_cols_per_proc)>=0) && ((coo-prank*num_cols_per_proc)<num_cols_per_proc)){
                vector <int> var;
                var.push_back(roo);
                var.push_back(coo);
                var.push_back(ino);
                var.push_back(flo);
                list_of_cols[coo-prank*num_cols_per_proc].push_back(var);
            }
            else{
                cout<<" Bad, in buffer to cols"<<endl;
            }
        }
        free(sendbfr);
        free(recvbfr);



        //sort cols
        for (i=0;i<list_of_cols.size();i++){
            // argsort one col according to float
            vector <pair<float, pair<int,int> > > tosort;
            for(j=0;j<list_of_cols[i].size();j++){
                tosort.push_back(make_pair(list_of_cols[i][j][2] ,make_pair(list_of_cols[i][j][0],j)));
            }
            sort(tosort.begin(),tosort.end());
            vector <vector<int> > newcol;
            for(j=0;j<list_of_cols[i].size();j++){
                vector<int> newel;
                newel.push_back(row_indices_in_cols[i][j]);
                newel.push_back(list_of_cols[i][tosort[j].second.second][1]);
                newel.push_back(list_of_cols[i][tosort[j].second.second][2]);
                newel.push_back(list_of_cols[i][tosort[j].second.second][3]);
                newcol.push_back(newel);
            }
            for(j=0;j<list_of_cols[i].size();j++){
                for(k=0;k<4;k++){
                    list_of_cols[i][j][k] = newcol[j][k];
                }
            }
            newcol = vector <vector<int> >();
            tosort =  vector <pair<float, pair<int,int> > >();
        }

        sendbfr = (int *) calloc(sum_recv_counts,sizeof(int));
        recvbfr = (int *) calloc(sum_send_counts,sizeof(int));

        added_to_s_bfr  = 0;
        //denotes how far I have travelled along different rows
        vector<int> col_trav(num_cols_per_proc,0);
        for(i=0;i<psize;i++){
            for(j=0;j<list_of_cols.size();j++){
                for(k=col_trav[j];k<list_of_cols[j].size();k++){
                    if((list_of_cols[j][k][0]>=(i*num_rows_per_proc)) && (list_of_cols[j][k][0]<(i+1)*num_rows_per_proc)){
                        sendbfr[added_to_s_bfr] = list_of_cols[j][k][0];
                        sendbfr[added_to_s_bfr+1] = list_of_cols[j][k][1];
                        sendbfr[added_to_s_bfr+2] = list_of_cols[j][k][2];
                        sendbfr[added_to_s_bfr+3] = list_of_cols[j][k][3];
                        added_to_s_bfr += 4;
                    }
                    else{
                        col_trav[j] = k;
                        break;
                    }
                }
            }
        }
        col_trav = vector<int>();

        // for(i=0;i<psize;i++){
        //     for(j=0;j<list_of_cols.size();j++){
        //         for(k=0;k<list_of_cols[j].size();k++){
        //             if((list_of_cols[j][k][0]>=(i*num_rows_per_proc)) && (list_of_cols[j][k][0]<(i+1)*num_rows_per_proc)){
        //                 sendbfr[added_to_s_bfr] = list_of_cols[j][k][0];
        //                 sendbfr[added_to_s_bfr+1] = list_of_cols[j][k][1];
        //                 sendbfr[added_to_s_bfr+2] = list_of_cols[j][k][2];
        //                 sendbfr[added_to_s_bfr+3] = list_of_cols[j][k][3];
        //                 added_to_s_bfr += 4;
        //             }
        //         }
        //     }
        // }

        
        MPI_Barrier(MPI_COMM_WORLD);
        // note that names here means opposite of the argument name (except buffer names)
        MPI_Alltoallv(sendbfr, recvcounts, rdispls, MPI_INT, recvbfr, sendcounts, sdispls, MPI_INT, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        
        list_of_rows = vector<vector< vector<int> > > ();
        list_of_rows.resize(num_rows_per_proc) ;
        list_of_cols = vector<vector< vector<int> > > ();

        // buffer to rows
        // sumsend counts is actually sum recv counts here
        for(i=0;i<sum_send_counts/4;i++){
            int ro = recvbfr[i*4];
            int co = recvbfr[i*4+1];
            int in = recvbfr[i*4+2];
            int fl = recvbfr[i*4+3];
            if(((ro-prank*num_rows_per_proc)>=0) && ((ro-prank*num_rows_per_proc)<num_rows_per_proc)){
                vector <int> var;
                var.push_back(ro);
                var.push_back(co);
                var.push_back(in);
                var.push_back(fl);
                list_of_rows[ro-prank*num_rows_per_proc].push_back(var);
            }
            else{
                cout<<" Bad, in buffer to rows"<<endl;
            }
        }
        free(sendbfr);
        free(recvbfr);
    }

    // cout<<"AFTER"<<endl;
    MPI_Barrier(MPI_COMM_WORLD);

    int trank = 0;
    FILE *of;
    while (trank < psize) {
        if (prank == trank) {
                // cout<<"Proc - "<<prank<<" Printing  My Rows"<<endl;
                // for (i=0;i<list_of_rows.size();i++){
                //     for(j=0;j<list_of_rows[i].size();j++){
                //         cout<<"row - "<<list_of_rows[i][j][0]<<" col - "<<list_of_rows[i][j][1]<< " i - "<<list_of_rows[i][j][2]<< " f - "<<*reinterpret_cast<float *>(&list_of_rows[i][j][3])<<endl;
                //     }
                //     cout<<endl<<"------------------------------------------------"<<endl;
                // }

            if (prank == 0){
                // ofstream clear_f(out_filename, ios::out | ios::trunc);
                of = fopen(out_filename,"wb");
                fflush(of);                
                fclose(of);
                // clear_f.close();
            }
            // ofstream my_output(out_filename,  ios::out | ios::binary | ios::app); 
            of = fopen(out_filename,"ab");            
            for (i=0;i<list_of_rows.size();i++){
                if(list_of_rows[i].size()>=1){
                    for(j=0;j<list_of_rows[i].size();j++){
                        float f = *reinterpret_cast<float *>( &(list_of_rows[i][j][3]) );
                        if(of!=NULL){
                            fwrite(&list_of_rows[i][j][0],sizeof(int),1,of);
                            fflush(of);
                            fwrite(&list_of_rows[i][j][1],sizeof(int),1,of);
                            fflush(of);                    
                            fwrite(&list_of_rows[i][j][2],sizeof(int),1,of);
                            fflush(of);                    
                            fwrite(&f,sizeof(float),1,of);
                            fflush(of);       
                        }      
                        else{
                            cout<<"File Pointer is null"<<endl;
                        }

                        // cout<<"Writing row - "<<list_of_rows[i][j][0]<<" col - "<<list_of_rows[i][j][1]<< " i - "<<list_of_rows[i][j][2]<< " f - "<<f<<endl;
                        // my_output << list_of_rows[i][j][0];
                        // my_output << list_of_rows[i][j][1];
                        // my_output << list_of_rows[i][j][2];
                        // my_output << f;
                        // my_output.write(list_of_rows[i][j][0], sizeof(int));
                        // my_output.write(list_of_rows[i][j][1], sizeof(int));
                        // my_output.write(list_of_rows[i][j][2], sizeof(int));                    
                        // my_output.write(f,sizeof(float));
                        // my_output.write(myrows[i][j].second,5);
                    }
                }
            }
            fclose(of);            
            // my_output.close();
        }
        trank ++;
        MPI_Barrier (MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}