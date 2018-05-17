
      
      i m_num_rows
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
                    val[j][k] = crecvbfr[ crdispls[proc_rank] + (maxlengthstring+1)*(i*column_counts[proc_rank] + col_to_proc_off[j] ) + k ];                    
                    // if(proc_rank != 0){
                    // val[j][k] = crecvbfr[ crdispls[proc_rank-1] + (maxlengthstring+1)*(i*column_counts[proc_rank] + col_to_proc_off[j] ) + k ];
                    // }
                    // else{
                    // val[j][k] = crecvbfr[ (maxlengthstring+1)*(i*column_counts[proc_rank] + col_to_proc_off[j] ) + k ];
                    // }
                }
                pair<float,char *> kv(inp_f,val[j]);
                row.push_back(kv);