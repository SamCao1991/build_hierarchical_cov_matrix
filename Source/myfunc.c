#include <stdio.h>
#include <time.h>
#include <math.h>

#include "aca.h"
#include "amatrix.h"
#include "cluster.h"
#include "clustergeometry.h"
#include "hmatrix.h"
#include "harith.h"

#define COLS 2 // which is the same as dim parameter

pamatrix expcovgen(const double geom[][COLS], uint dim, pccluster clus1, pccluster clus2, const double beta);

phmatrix myBuildHmatrix_exp_0822(const double geom[][COLS], uint dim, pccluster root,const double beta,const double tol)
{
    phmatrix hm;
    if(root->sons == 0)
    {
        hm = new_hmatrix(root,root);
        hm->f = expcovgen(geom, dim, root, root, beta);
    }else
    {
        assert(root->sons == 2);
        hm = new_super_hmatrix(root,root,2,2);
        ref_hmatrix(hm->son,myBuildHmatrix_exp_0822(geom, dim, root->son[0],beta,tol));
        ref_hmatrix(hm->son+1,new_rk_hmatrix(root->son[1],root->son[0],1));
        pamatrix fullmt_temp = expcovgen(geom, dim, root->son[1],root->son[0],beta);
        decomp_fullaca_rkmatrix(fullmt_temp,tol,NULL,NULL,hm->son[1]->r);
        del_amatrix(fullmt_temp);
        ref_hmatrix(hm->son+2,new_rk_hmatrix(root->son[0],root->son[1],1));
        fullmt_temp = expcovgen(geom, dim, root->son[0],root->son[1],beta);
        decomp_fullaca_rkmatrix(fullmt_temp,tol,NULL,NULL,hm->son[2]->r);
        del_amatrix(fullmt_temp);
        ref_hmatrix(hm->son+3,myBuildHmatrix_exp_0822(geom, dim, root->son[1],beta,tol));
    }
    update_hmatrix(hm);
    return hm;
}

void myOutputB(pchmatrix hm, FILE *fp)
{
    if(hm->f != NULL)
    {
        int i,j;
        for(i = 0; i < hm->f->rows; i++)
        {
            for(j = 0; j < hm->f->cols; j++)
            {
                fprintf(fp,"%lf ",hm->f->a[i+j*hm->f->ld]);
            }
            fprintf(fp,"\n");
        }
    }
    else
    {
        myOutputB(hm->son[0],fp);
        myOutputB(hm->son[3],fp);
    }
}

void myOutputUV(pchmatrix hm, FILE *fp, int q)
{
    if(hm->f == NULL)
    {
        int i,j;
        myOutputUV(hm->son[0],fp,q+1);
        fprintf(fp,"%d\n%d\n%d\n%d\n%d\n",hm->son[1]->rc->idx[0],hm->son[1]->cc->idx[0],hm->son[1]->rc->size,q,hm->son[1]->r->k);
        for(i = 0; i < hm->son[1]->rc->size; i++)
        {
            for(j = 0; j < hm->son[1]->r->k; j++)
            {
                fprintf(fp,"%lf ",hm->son[1]->r->A.a[i+j*hm->son[1]->r->A.ld]);
            }
            fprintf(fp,"\n");
        }
        for(i = 0; i < hm->son[1]->cc->size; i++)
        {
            for(j = 0; j < hm->son[1]->r->k; j++)
            {
                fprintf(fp,"%lf ",hm->son[1]->r->B.a[i+j*hm->son[1]->r->B.ld]);
            }
            fprintf(fp,"\n");
        }
        myOutputUV(hm->son[3],fp,q+1);
    }
}

pamatrix expcovgen(const double geom[][COLS], uint dim, pccluster clus1, pccluster clus2, const double beta)
{
    pamatrix fullmt = new_amatrix(clus1->size,clus2->size);
    init_amatrix(fullmt,clus1->size,clus2->size);
    int i,j,k;
    double dist;
    for(i=0;i<clus1->size;i++)
    {
        for(j=0;j<clus2->size;j++)
        {
            dist = 0;
            for(k=0;k<dim;k++)
                dist += (geom[clus1->idx[i]][k] - geom[clus2->idx[j]][k]) * (geom[clus1->idx[i]][k] - geom[clus2->idx[j]][k]);
            dist = sqrt(dist);
            fullmt->a[i+fullmt->ld*j] = exp(- dist / beta);
        }
    }
    return fullmt;
}

void output_exp_hie_cov(const double geom[][COLS], uint N, uint dim, uint bsz, double beta, double tol)
{
    pclustergeometry clustergeom;
    pcluster clus;
    phmatrix hm;
    clustergeom = new_clustergeometry(1,N);
    uint *idx = (uint *)malloc(sizeof(int) * N);
    int i = 0;
    while(i < N)
    {
        clustergeom->x[i][0] = i;
        idx[i] = i;
        i++;
    }
    update_point_bbox_clustergeometry(clustergeom,N,idx);
    clus = build_simsub_cluster(clustergeom,N,idx,bsz);
    hm = myBuildHmatrix_exp_0822(geom,dim,clus,beta,tol);
    choldecomp_hmatrix(hm,0,tol);
    FILE *fp;
    fp = fopen("B.txt","w+");
    fprintf(fp,"%d\n",bsz);
    myOutputB(hm,fp);
    fclose(fp);
    fp = fopen("UV.txt","w+");
    myOutputUV(hm,fp,0);
    fclose(fp);

    cairo_t *ppdf;
    ppdf = new_cairopdf("hmatrix.pdf",1024.0,1024.0);
    draw_cairo_hmatrix(ppdf,hm,false,0);
    cairo_destroy(ppdf);

    free(idx);
    del_hmatrix(hm);
    del_cluster(clus);
    del_clustergeometry(clustergeom);
}


