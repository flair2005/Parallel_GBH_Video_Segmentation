/*
Original Code From:
Copyright (C) 2006 Pedro Felzenszwalb
Modifications (may have been made) Copyright (C) 2011,2012 
  Chenliang Xu, Jason Corso.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

/* Implements node merging criteria. */

#ifndef SEGMENT_GRAPH_MULTI
#define SEGMENT_GRAPH_MULTI

#include <algorithm>
#include <cmath>
#include <vector>
#include <omp.h>
#include "pnmfile.h"
#include "disjoint-set.h"
#include "disjoint-set-s.h"
#include "segment-graph-s.h"
#include "edges.h"

using namespace std;

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

// random color
rgb random_rgb() {
	rgb c;
 	c.r = (uchar) random();
        c.g = (uchar) random();
        c.b = (uchar) random();
        return c;
}
/*
bool operator<(const edge &a, const edge &b) {

	return a.w < b.w;

}
*/
// dissimilarity measure between pixels
static inline float diff(image<float> *r, image<float> *g, image<float> *b,
                         int x1, int y1, int x2, int y2) {
  return sqrt(square(imRef(r, x1, y1)-imRef(r, x2, y2)) +
              square(imRef(g, x1, y1)-imRef(g, x2, y2)) +
              square(imRef(b, x1, y1)-imRef(b, x2, y2)));
}

// process every image with graph-based segmentation
void gb(universe *mess, image<float> *smooth_r, image<float> *smooth_g, image<float> *smooth_b,
        int width, int height, edge *edges, float c, int s_index, int e_index, int level, 
        vector<edge>* edges_remain, int num) {
  //int num = 0;
  //  edge_s *edges = new edge_s[width*height*4];
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      if (x < width-1) {
        edges[num].a = y * width + x;
        edges[num].b = y * width + (x+1);
        edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y);
        num++;
      }
      if (y < height-1) {
        edges[num].a = y * width + x;
        edges[num].b = (y+1) * width + x;
        edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x, y+1);
        num++;
      }
      if ((x < width-1) && (y < height-1)) {
        edges[num].a = y * width + x;
        edges[num].b = (y+1) * width + (x+1);
        edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y+1);
        num++;
      }
      if ((x < width-1) && (y > 0)) {
        edges[num].a = y * width + x;
        edges[num].b = (y-1) * width + (x+1);
        edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y-1);
        num++;
      }
    }
  }

  universe_s *u = segment_graph_s(width*height, num, edges, c, edges_remain);
 
  for (int i = s_index; i < e_index; ++i) 
    mess->set_in_level(i, level, u->find(i-s_index), u->rank(i-s_index), u->size(i-s_index), u->mst(i-s_index)); 
}


/* pixel level minimum spanning tree merge */
void segment_graph(universe *mess, vector<edge>* edges_remain, edge *edges, 
		int num_edges, float c, int width, int height, int level,
                image<float> *smooth_r[], image<float> *smooth_g[], image<float> *smooth_b[]) {
	// new vector containing remain edges
	edges_remain->clear();

	printf("Start segmenting graph in parallel.\n");
	int th_id;
	vector<edge>* edges_remain0 = new vector<edge>();
	vector<edge>* edges_remain1 = new vector<edge>();
	vector<edge>* edges_remain2 = new vector<edge>();
	vector<edge>* edges_remain3 = new vector<edge>();
	vector<edge>* edges_remain4 = new vector<edge>();
	vector<edge>* edges_remain5 = new vector<edge>();
	vector<edge>* edges_remain6 = new vector<edge>();
	vector<edge>* edges_remain7 = new vector<edge>();
//	int sfpxl_num = width * height; // single frame pixel number
        #pragma omp parallel private(th_id) 
	{
  	  th_id = omp_get_thread_num();
          switch(th_id) {
            case 0: 
            {
	      int num0 = 0;
              edge *edges0 = new edge[width*height*4];
              gb(mess, smooth_r[0], smooth_g[0], smooth_b[0], width, height, edges0, c, 0, width*height, level, edges_remain0, num0);
            }
	    break;
            case 1: 
            {
	      int num1 = 0;
              edge *edges1 = new edge[width*height*4];
              gb(mess, smooth_r[1], smooth_g[1], smooth_b[1], width, height, edges1, c, width*height, 2*width*height, level, edges_remain1, num1);            
            }
            break;
            case 2: 
       	    {
	      int num2 = 0;
 	      edge *edges2 = new edge[width*height*4];
	      gb(mess, smooth_r[2], smooth_g[2], smooth_b[2], width, height, edges2, c, 2*width*height, 3*width*height, level, edges_remain2, num2);            
            }
            break;
            case 3: 
            {
	      int num3 = 0;
	      edge *edges3 = new edge[width*height*4];
	      gb(mess, smooth_r[3], smooth_g[3], smooth_b[3], width, height, edges3, c, 3*width*height, 4*width*height, level, edges_remain3, num3);            
            }
            break;
            case 4: 
            {
	      int num4 = 0;
	      edge *edges4 = new edge[width*height*4];
	      gb(mess, smooth_r[4], smooth_g[4], smooth_b[4], width, height, edges4, c, 4*width*height, 5*width*height, level, edges_remain4, num4);            
            }
      	    break;
            case 5: 
            {
	      int num5 = 0;
	      edge *edges5 = new edge[width*height*4];
              gb(mess, smooth_r[5], smooth_g[5], smooth_b[5], width, height, edges5, c, 5*width*height, 6*width*height, level, edges_remain5, num5);            
	    }
            break;
            case 6: 
            {
	      int num6 = 0;
	      edge *edges6 = new edge[width*height*4];
	      gb(mess, smooth_r[6], smooth_g[6], smooth_b[6], width, height, edges6, c, 6*width*height, 7*width*height, level, edges_remain6, num6);            
       	    }
       	    break;
      	    case 7: 
       	    {
	      int num7 = 0;
	      edge *edges7 = new edge[width*height*4];
	      gb(mess, smooth_r[7], smooth_g[7], smooth_b[7], width, height, edges7, c, 7*width*height, 8*width*height, level, edges_remain7, num7);            
       	    }
       	    break;
     
       	    default: break;
	  }
    	}
//        for (int i = 0; i < num0; ++i) {
//          edge *pedge = &edges0[i];
//	}
	vector<edge>::iterator it;
        for ( it = edges_remain0->begin() ; it < edges_remain0->end(); it++ )
          edges_remain->push_back(*it); 
        for ( it = edges_remain1->begin() ; it < edges_remain1->end(); it++ )
          edges_remain->push_back(*it); 
        for ( it = edges_remain2->begin() ; it < edges_remain2->end(); it++ )
          edges_remain->push_back(*it); 
        for ( it = edges_remain3->begin() ; it < edges_remain3->end(); it++ )
          edges_remain->push_back(*it); 
        for ( it = edges_remain4->begin() ; it < edges_remain4->end(); it++ )
          edges_remain->push_back(*it); 
        for ( it = edges_remain5->begin() ; it < edges_remain5->end(); it++ )
          edges_remain->push_back(*it); 
        for ( it = edges_remain6->begin() ; it < edges_remain6->end(); it++ )
          edges_remain->push_back(*it); 
        for ( it = edges_remain7->begin() ; it < edges_remain7->end(); it++ )
          edges_remain->push_back(*it); 
	sort(edges_remain->begin(), edges_remain->end());
}
	
/* region graph level minimum spanning tree merge */
void segment_graph_region(universe *mess, vector<edge>* edges_remain, 
		vector<edge>* edges_region, float c_reg, int level) {  
	edges_remain->clear();
	sort(edges_region->begin(), edges_region->end());
//	printf("Finished edge region sorting.\n");
	printf("Edge region size is %d.\n", (int)edges_region->size());
	for (int i=0; i < (int) edges_region->size(); i++) {
//		printf("Operation+1\n");
//		float max_w = 0;
//	        if (max_w < edges_region->at(i).w)
//	          max_w = edges_region->at(i).w;
//	        printf("Maximum region edge weight is %f.\n", max_w); 
		int a = mess->find_in_level(edges_region->at(i).a, level);
		int b = mess->find_in_level(edges_region->at(i).b, level);
//                float a_mst = mess->get_mst(a);
//                float b_mst = mess->get_mst(b);
//                printf("Maximum mst is %f.\n",(a_mst>b_mst?a_mst:b_mst));
		if (a != b) {
			if ((edges_region->at(i).w
					<= mess->get_mst(a) + (c_reg / mess->get_size(a)))
					&& (edges_region->at(i).w
							<= mess->get_mst(b) + (c_reg / mess->get_size(b)))) {
				if (mess->join(a, b, edges_region->at(i).w, level) == 1)
					edges_remain->push_back(edges_region->at(i));
			} else {
				edges_remain->push_back(edges_region->at(i));
			}
//		printf("Operation+1\n");
		}

	}
}
#endif /* SEGMENT_GRAPH_H */
