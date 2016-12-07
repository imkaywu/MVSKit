// input data:
//   - 2-view images
//   - images under different illuminations
//   - normal of 2 views
//   - initial discrete dense correspondence
//   - camera pamameters (intrinsic and extrinsic)

// output data:
//   - continuous dense correspondence (no need for rectification)
//   - 3d point cloud

/* -------------------------------------------------------------------------
 BITMAP: Minimal image class
 ------------------------------------------------------------------------- */
#include "pm_mvps.hpp"

int num_iter = 5; // 2 for each direction
int patch_hw = 3, patch_w = 2 * patch_hw + 1;
int cutoff = INT_MAX;

void check_im()
{
    if (system("identify > null.txt") != 0) {
        fprintf(stderr, "ImageMagick must be installed, and 'convert' and 'identify' must be in the path\n"); exit(1);
    }
}

BITMAP *load_bitmap(const char *filename)
{
    check_im();
    char rawname[256], txtname[256];
    strcpy(rawname, filename);
    strcpy(txtname, filename);
    if (!strstr(rawname, ".")) { fprintf(stderr, "Error reading image '%s': no extension found\n", filename); exit(1); }
    sprintf(strstr(rawname, "."), ".raw");
    sprintf(strstr(txtname, "."), ".txt");
    char buf[256];
    sprintf(buf, "convert %s rgba:%s", filename, rawname);
    if (system(buf) != 0) { fprintf(stderr, "Error reading image '%s': ImageMagick convert gave an error\n", filename); exit(1); }
    sprintf(buf, "identify -format \"%%w %%h\" %s > %s", filename, txtname);
    if (system(buf) != 0) { fprintf(stderr, "Error reading image '%s': ImageMagick identify gave an error\n", filename); exit(1); }
    FILE *f = fopen(txtname, "rt");
    if (!f) { fprintf(stderr, "Error reading image '%s': could not read output of ImageMagick identify\n", filename); exit(1); }
    int w = 0, h = 0;
    if (fscanf(f, "%d %d", &w, &h) != 2) { fprintf(stderr, "Error reading image '%s': could not get size from ImageMagick identify\n", filename); exit(1); }
    fclose(f);
    f = fopen(rawname, "rb");
    BITMAP *ans = new BITMAP(w, h);
    unsigned char *p = (unsigned char *) ans->data;
    for (int i = 0; i < w*h*4; i++) {
        int ch = fgetc(f);
        if (ch == EOF) { fprintf(stderr, "Error reading image '%s': raw file is smaller than expected size %dx%dx4\n", filename, w, h, 4); exit(1); }
        *p++ = ch;
    }
    fclose(f);
    return ans;
}

void save_bitmap(BITMAP *bmp, const char *filename)
{
    check_im();
    char rawname[256];
    strcpy(rawname, filename);
    if (!strstr(rawname, ".")) { fprintf(stderr, "Error writing image '%s': no extension found\n", filename); exit(1); }
    sprintf(strstr(rawname, "."), ".raw");
    char buf[256];
    FILE *f = fopen(rawname, "wb");
    if (!f) { fprintf(stderr, "Error writing image '%s': could not open raw temporary file\n", filename); exit(1); }
    unsigned char *p = (unsigned char *) bmp->data;
    for (int i = 0; i < bmp->w*bmp->h*4; i++) {
        fputc(*p++, f);
    }
    fclose(f);
    sprintf(buf, "convert -size %dx%d -depth 8 rgba:%s %s", bmp->w, bmp->h, rawname, filename);
    if (system(buf) != 0) { fprintf(stderr, "Error writing image '%s': ImageMagick convert gave an error\n", filename); exit(1); }
}

double *load_txt(const char* fname) // used to load matrix from text files.
{
    return 0;
}

/* -------------------------------------------------------------------------
 PatchMatch Stereo
 ------------------------------------------------------------------------- */
// matchl2r, matchr2l: W x H x 2 matrix, 1st channel stores x-coordinate, 2nd channel stores y-coordinate
// nl, nr: W x H x 3 matrix
// Pl, Pr: projection matrix for left and right cameras
void patchmatch_stereo_main(char *image_base, char *mask_base, int num_img, double *matchl2r, double *matchr2l, double *nl, double *nr, vector<Matrix3Xd> &P, double *depth)
{
    vector<BITMAP *> a(num_img), b(num_img); // don't forget to delete
    vector<BITMAP *> mask_ref(2), mask_tar(2);
    std::string image_name, mask_name;
    for(int i = 0; i < num_img; ++i)
    {
        if(i < 10)
        {
            image_name = string(image_base) + "_0_0" + std::to_string(i) + ".BMP";
            a[i] = load_bitmap(image_name.c_str());
            image_name = string(image_base) + "_1_0" + std::to_string(i) + ".BMP";
            b[i] = load_bitmap(image_name.c_str());
        }
        else
        {
            image_name = string(image_base) + "_0_" + std::to_string(i) + ".BMP";
            a[i] = load_bitmap(image_name.c_str());
            image_name = string(image_base) + "_1_" + std::to_string(i) + ".BMP";
            b[i] = load_bitmap(image_name.c_str());
        }
    }
    mask_name = string(mask_base) + "_ref_0.BMP";
    mask_ref[0] = load_bitmap(mask_name.c_str());
    mask_name = string(mask_base) + "_ref_1.BMP";
    mask_ref[1] = load_bitmap(mask_name.c_str());
    mask_name = string(mask_base) + "_target_0.BMP";
    mask_tar[0] = load_bitmap(mask_name.c_str());
    mask_name = string(mask_base) + "_target_1.BMP";
    mask_tar[1] = load_bitmap(mask_name.c_str());
    
    patchmatch_stereo(a, b, mask_ref, mask_tar, matchl2r, matchr2l, nl, nr, P, depth);
}

void patchmatch_stereo(vector<BITMAP *> &a, vector<BITMAP *> &b, vector<BITMAP *> &mask_ref, vector<BITMAP *> &mask_tar, double *matchl2r, double *matchr2l, double *nl, double *nr, vector<Matrix3Xd> &P, double *depth)
{
    int aw = a[0]->w, ah = a[0]->h;
    int bw = b[0]->w, bh = b[0]->h;
    int aew = aw - patch_w + 1, aeh = ah - patch_w + 1;
    int bew = bw - patch_w + 1, beh = bh - patch_w + 1;
    int reso = aw * ah;
    
    double bx, by;
    double *dis = new double[aw * ah]; // dissimilarity score
    memset(dis, 0, sizeof(double) * reso);
    for (int ay = 0; ay < ah; ++ay)
    {
        for (int ax = 0; ax < aw; ++ax)
        {
            if(*mask_tar[0][ay][ax] > 0)
            {
                bx = matchl2r[ay * aw + ax + 0 * reso];
                by = matchl2r[ay * aw + ax + 1 * reso];
                dis[ay * aw + ax] = dist(a, b, nl, nr, ax, ay, bx, by);
            }
        }
    }
    
    for (int iter = 0; iter < num_iter; ++iter)
    {
        int ystart = patch_hw - 1, yend = ah - patch_hw + 1, ychange = 1;
        int xstart = patch_hw - 1, xend = aw - patch_hw + 1, xchange = 1;
        if(iter % 2 == 1)
        {
            ystart = yend - patch_hw + 1; yend = -1; ychange = -1;
            xstart = xend - patch_hw + 1; xend = -1; xchange = -1;
        }
        
        for (int ay = ystart; ay != yend; ay += ychange)
        {
            for (int ax = xstart; ax != xend; ax += xchange)
            {
                if(*mask_tar[0][ay][ax] <= 0)
                    continue;
                
                // current guess
                double xbest = matchl2r[ay * aw + ax + 0 * reso]; // continuous correspondence
                double ybest = matchl2r[ay * aw + ax + 1 * reso];
                double dbest = dis[ay * aw + ax];
                
                // spatial propagation
                if((unsigned) (ax - xchange) < (unsigned) aew && *mask_tar[0][ay][ax - xchange] > 0)
                {
                    int xp = ax - xchange, yp = ay; // left/right point
                    improve_guess(a, b, nl, nr, ax, ay, xp, yp, xbest, ybest, dbest);
                }
                
                if((unsigned) (ay - ychange) < (unsigned) aeh && *mask_tar[0][ay - ychange][ax] > 0)
                {
                    int xp = ax, yp = ay - ychange; // top/bottom point
                    improve_guess(a, b, nl, nr, ax, ay, xp, yp, xbest, ybest, dbest);
                }
                
                // view propagation
                for (int dy = -patch_hw; dy < patch_hw; ++dy)
                {
                    for (int dx = -patch_hw; dx < patch_hw; ++dx)
                    {
                        bx = matchr2l[int(ybest + dy + 0.5) * bw + int(xbest + dx + 0.5) + 0 * reso];
                        by = matchr2l[int(ybest + dy + 0.5) * bw + int(xbest + dx + 0.5) + 1 * reso];
                        if(int(bx + 0.5) == ax && int(by + 0.5) == ay) // round to the closest integer
                        {
                            improve_guess(a, b, nl, nr, ax, ay, xbest + dx, ybest + dy, xbest, ybest, dbest);
                        }
                    }
                }
                
                // plane refinement
                Vector3d X;
                Vector2d u, v;
                u << ax, ay;
                v << xbest, ybest;
                vector<Vector2d> x(2);
                x[0] = u, x[1] = v;
                vector<Vector2d> imsize(2);
                imsize[0](0) = ah, imsize[0][1] = aw;
                imsize[1](0) = bh, imsize[1][1] = bw;
                X = triangulate_nonlin(P, imsize, x);
                double d = pt2depth(P[0], X);
                double thre = 1e-3, shrink = 0.5;
                double dd = 0.5 * d, drand, dtmp;
                while(1)
                {
                    if(dd < thre * d)
                        break;
                    drand = uniform_random(dd); // a random number within (0, dd)
                    
                    dtmp = d + drand;
                    X = depth2pt(P[0], u, dtmp);
                    v = reproj(P[1], X);
                    improve_guess(a, b, nl, nr, ax, ay, bx, by, xbest, ybest, dbest);
                    
                    dtmp = d - drand;
                    X = depth2pt(P[0], u, dtmp);
                    v = reproj(P[1], X);
                    improve_guess(a, b, nl, nr, ax, ay, bx, by, xbest, ybest, dbest);
                    
                    dd = shrink * dd;
                }
            }
        }
    }
    delete dis;
}

double dist(vector<BITMAP *> &a, vector<BITMAP *> &b, double *nl, double *nr, double ax, double ay, double bx, double by)
{
    double ans = 0;
    int w = a[0]->w, h = a[0]->h, reso = w * h;
    ans = ans + nl[int(ay + 0.5) * w + int(ax + 0.5) + 0 * reso] * nr[int(by * w + 0.5) + int(bx + 0.5) + 0 * reso]; // normal cost
    ans = ans + nl[int(ay + 0.5) * w + int(ax + 0.5) + 1 * reso] * nr[int(by * w + 0.5) + int(bx + 0.5) + 1 * reso];
    ans = ans + nl[int(ay + 0.5) * w + int(ax + 0.5) + 2 * reso] * nr[int(by * w + 0.5) + int(bx + 0.5) + 2 * reso];
    
    int num_img = static_cast<int>(a.size());
    for (int n = 0; n < num_img; ++n)
    {
        for (int dy = 0; dy < patch_w; ++dy)
        {
            int *arow = &(*a[n])[int(ay + 0.5)+dy][int(ax + 0.5)];
            int *brow = &(*b[n])[int(by + 0.5)+dy][int(bx + 0.5)];
            for (int dx = 0; dx < patch_w; ++dx)
            {
                int ac = arow[dx];
                int bc = brow[dx];
                int dr = (ac&255)-(bc&255);
                int dg = ((ac>>8)&255)-((bc>>8)&255);
                int db = (ac>>16)-(bc>>16);
                ans += dr*dr + dg*dg + db*db;
            }
            if (ans >= cutoff) { return cutoff; }
        }
    }
    return ans;
}

void improve_guess(vector<BITMAP *> &a, vector<BITMAP *> &b, double *nl, double *nr, double ax, double ay, double bx, double by, double xbest, double ybest, double dbest)
{
    double d = dist(a, b, nl, nr, ax, ay, bx, by);
    if(d < dbest)
    {
        xbest = bx;
        ybest = by;
        dbest = d;
    }
}

double uniform_random(double range)
{
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, range);
    return distribution(generator);
}

void compt_grad(BITMAP &mask_tar, double *depth) // for normal cost function
{
    
}


