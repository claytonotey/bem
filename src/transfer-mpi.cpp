#include "mpi.h"
#include "HeatTransfer.h"

using namespace std;
#include <queue>

#define TAG_JOB 1
#define TAG_RET 2

class FJob {
public:
  Real buf[8];
  FJob(Real a, const Vector3 &gap, Real feV, Real T1, Real T2) {
    buf[0] = a;
    buf[1] = gap.x;
    buf[2] = gap.y;
    buf[3] = gap.z;
    buf[4] = feV;
    buf[5] = T1;
    buf[6] = T2;
  }    
};


class Servant {
public:
  RWGTree *top;
  RWGTree *s;
  RWGTree *t;
  int masterProc;
  Options &options;

  Servant(int masterProc, Options &options) : 
    masterProc(masterProc), options(options) {}
  
  void start() {
    Real buf[7];
    MPI_Status status;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   
    while(true) {

      MPI_Recv(buf,5,MPI_DOUBLE,masterProc,TAG_JOB,MPI_COMM_WORLD,&status);

      if(buf[0] < 0) break;

      Real a = buf[0];
      Vector3 gap(buf[1],buf[2],buf[3]);
      Real eV = buf[4];
      Real T1 = buf[5];
      Real T2 = buf[6];
      char body1x3d[32] = "body1.x3d";
      char body2x3d[32] = "body2.x3d";

      parseShapes(top, s, t, body1x3d, body2x3d, gap, options);
      Real S = CalcHeatTransfer(a, eV, T1, T2, top, s, t, options);
      Real out[1];
      out[0] = S;
      MPI_Send(out,1,MPI_DOUBLE,masterProc,TAG_RET,MPI_COMM_WORLD);    
    }
    fprintf(stdout,"servant %d exiting...\n",rank);
  }
  
};


class Master {
public:
  vector<Real> &as;
  vector<Vector3> &gaps;
  vector<Point2> &temps;
  vector<Real> &freqs;
  Array<FJob*> fjobs;
  Array<int> procJ;
  Array<Real> procResult;
  Array<MPI_Request> procRequest;
  queue<int> procQ;
  queue<int> fjobQ;
  int bufsize;
  int nprocs;
  int na;
  int nf;
  int nfreqs;
  int ngaps;
  int nfjobs;
  int nrequests;
  int ntemps;
  Master(vector<Real> &as, vector<Vector3> &gaps, vector<Point2> &temps, vector<Real> &freqs) : as(as), gaps(gaps), temps(temps), freqs(freqs) {}

  void start()
  {
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Status status;

    na = as.size();
    nf = freqs.size();
    ngaps = gaps.size();
    ntemps = temps.size();

    nfjobs = nf*ngaps*ntemps;
    fjobs.resize(nfjobs);
    procRequest.resize(nprocs);
    bufsize = 1;
    procResult.resize(nprocs*bufsize);
    procJ.resize(nprocs);
    nrequests = 0;

    FILE *out = fopen("output","w");

    int j = 0;
    for(int m=0; m<ntemps; m++) {
      Point2 T = temps[m];
      Real T1 = T.x;
      Real T2 = T.y;
      for(int l=0; l<na; l++) {
        Real a = as[l];
        for(int k=0; k<ngaps; k++) {
          Vector3 &gap = gaps[k];
          for(int i=0; i<nf; i++) {
            Real feV = freqs[i];
            FJob *fjob = new FJob(a,gap,feV,T1,T2);
            fjobQ.push(j);
            fjobs[j] = fjob;
            j++;
          }
        }
      }
    }
    
    for(int i=0; i<nprocs; i++) {
      procQ.push(i);
    }


    while(nfjobs) {
      while(!procQ.empty() && !fjobQ.empty()) {
        int proc = procQ.front(); procQ.pop();
        FJob *fjob;

        int j;
        do {
          j = fjobQ.front(); fjobQ.pop();
          fjob = fjobs[j];
        } while(fjob == NULL && !fjobQ.empty());

        if(fjob == NULL) break;

        MPI_Ssend(fjob->buf,5,MPI_DOUBLE,proc,TAG_JOB,MPI_COMM_WORLD);
        MPI_Irecv(procResult+proc*bufsize,bufsize,MPI_DOUBLE,proc,TAG_RET,MPI_COMM_WORLD,procRequest+proc);
        procJ[proc] = j;
        nrequests = min(nprocs,nrequests+1);
      }

      //fprintf(stderr,"%d requests\n",nrequests);
      int procDone;
      MPI_Waitany(nrequests,procRequest,&procDone,&status);

      if(procDone != MPI_UNDEFINED) {
        int j = procJ[procDone];
        FJob *fjob = fjobs[j];

        Real a = fjob->buf[0];
        Real feV = fjob->buf[1];
        Real gap = fjob->buf[2];
        Real T1 = fjob->buf[3];
        Real T2 = fjob->buf[4];
        Real S = procResult[procDone];
        fprintf(out,"%.15g %.15g %.15g %.15g %.15g %.14g\n",a,gap,T1,T2,feV,S);
        fflush(out);
        
        delete fjob;
        fjobs[j] = NULL;
        nfjobs--;
      }
      procQ.push(procDone);
    }
  
    for(int i=0; i<nprocs; i++) {
      Real buf[5];
      buf[0] = -1;
      MPI_Ssend(buf,5,MPI_DOUBLE,i,TAG_JOB,MPI_COMM_WORLD);
    }
    
    fprintf(stdout,"master %d exiting...\n",rank);
    fclose(out);
  }
};
  
void *masterThreadCB(void *data) {
  
  Master *master = (Master*)data;
  master->start();
  return NULL;
}


int main(int argc, char **argv)
{
  initReferenceCountingPointerLock();

  int rank;
  int provided;
  MPI_Init_thread(NULL,NULL,MPI_THREAD_MULTIPLE,&provided);

  int size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  vector<Real> as;
  vector<Vector3> gaps;
  vector<Real> freqs;
  vector<Point2> temps;


  char body1x3d[32] = "body1.x3d";
  char body2x3d[32] = "body2.x3d";

  pushFileLines("a.txt",as);
  pushFileLines("gaps.txt",gaps);
  pushFileLines("freqs.txt",freqs);
  pushFileLines("temps.txt",temps);

  Options options;
  if(rank == 0) {
    Master master(as,gaps,temps,freqs);
    pthread_t masterThread;
    pthread_create(&masterThread, NULL, masterThreadCB, &master);
    Servant servant(0, options);
    servant.start();
    pthread_join(masterThread,NULL);
  } else {
    Servant servant(0, options);
    servant.start();
  }

  MPI_Finalize();
  destroyReferenceCountingPointerLock();
}

