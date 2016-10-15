/*
 ============================================================================
 Name        : raycast.c
 Author      : Anthony Black
 Description : CS430 Project 2: Raycaster
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

//#define DEBUG 1 //uncomment to see print statements   
//STRUCTURES
// Plymorphism in C
typedef struct Object{
  int type; // 0 = sphere, 1 = plane
  double color[3];
  double position[3];
  union {
    struct {
      double radius;
    } sphere;
    struct {
      double normal[3];
    } plane;
  };
} Object;

typedef struct Camera{
  double width;
  double height;
} Camera;

typedef struct Pixel{
  unsigned char r, b, g;
}Pixel;

//PROTOTYPE DECLARATIONS

int next_c(FILE* json);

void expect_c(FILE* json, int d);

void skip_ws(FILE* json);

char* next_string(FILE* json);

double next_number(FILE* json);

double* next_vector(FILE* json);

void read_scene(char* filename, Camera* camera, Object** objects);

double sqr(double v);

void normalize(double* v);

double sphere_intersection(double* Ro, double* Rd, double* C, double r);

double plane_intersection(double* Ro, double* Rd, double* P, double* N);

void generate_scene(Camera* camera, Object** objects, Pixel* buffer, int width, int height);

void write_p3(Pixel *buffer, FILE *output_file, int width, int height, int max_color);

//Global variable for tracking during reading of JSON file, to report errors.
int line = 1;

//FUNCTIONS

int main(int argc, char *argv[]) {
  //ensures the correct number are passed in
  #ifdef DEBUG
    printf("Checking arguments...\n");
  #endif
  if (argc != 5){
    fprintf(stderr, "Error: Insufficient Arguments. Arguments provided: %d.\n", argc);
  }
  #ifdef DEBUG
    printf("Getting width and height...\n");
  #endif
  int width = atoi(argv[1]);
  int height = atoi(argv[2]);
  //check for positive width
  if (width <= 0){
    fprintf(stderr, "Error: Non-positive width provided.\n");
  }
  //check for positive height
  if (height <= 0){
    fprintf(stderr, "Error: Non-positive height provided.\n");
  }
  //create array of objects
  #ifdef DEBUG
    printf("Allocating memory...\n");
  #endif
  Object **objects;
  objects = malloc(sizeof(Object*)*128);
  //create camera object
  Camera *camera;
  camera = (Camera *)malloc(sizeof(Camera));
  //create buffer for image
  Pixel *buffer; 
  buffer = (Pixel *)malloc(width*height*sizeof(Pixel));
  #ifdef DEBUG
    printf("Reading scene...\n");
  #endif
  read_scene(argv[3], camera, objects);
  #ifdef DEBUG
    printf("Generating scene...\n");
  #endif
  generate_scene(camera, objects, buffer, width, height);
  #ifdef DEBUG
    printf("Opening output file...\n");
  #endif
  FILE* output_file = fopen(argv[4], "w");
  //error handling for failure to open output file
  if (output_file == NULL){
    fprintf(stderr, "Error: Unexpectedable to open output file.\n");
    exit(1);
  }
  #ifdef DEBUG
    printf("Creating image...\n");
  #endif
  write_p3(buffer, output_file, width, height, 255);
  return EXIT_SUCCESS;
}

int next_c(FILE* json) {
  // next_c() wraps the getc() function and provides error checking and line
  // number maintenance
  int c = fgetc(json);
  #ifdef DEBUG
    printf("next_c: '%c'\n", c);
  #endif
  if (c == '\n') {
    line += 1;
  }
  if (c == EOF) {
    fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", line);
    exit(1);
  }
  return c;
}


void expect_c(FILE* json, int d) {
  // expect_c() checks that the next character is d.  If it is not it emits
  // an error.
  int c = next_c(json);
  if (c == d) return;
  fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, line);
  exit(1);    
}


void skip_ws(FILE* json) {
  // skip_ws() skips white space in the file.
  int c = next_c(json);
  while (isspace(c)) {
    c = next_c(json);
  }
  ungetc(c, json);
}


char* next_string(FILE* json) {
  // next_string() gets the next string from the file handle and emits an error
  // if a string can not be obtained. 
  char buffer[129];
  int c = next_c(json);
  if (c != '"') {
    fprintf(stderr, "Error: Expected string on line %d.\n", line);
    exit(1);
  }  
  c = next_c(json);
  int i = 0;
  while (c != '"') {
    if (i >= 128) {
      fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
      exit(1);      
    }
    if (c == '\\') {
      fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
      exit(1);      
    }
    if (c < 32 || c > 126) {
      fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
      exit(1);
    }
    buffer[i] = c;
    i += 1;
    c = next_c(json);
  }
  buffer[i] = 0;
  return strdup(buffer);
}

double next_number(FILE* json) {
  double value;
  int count = fscanf(json, "%lf", &value);
  if (count != 1){
    fprintf(stderr, "Error: Failed to read number on line %d.\n", line);
    exit(1);
  }
  return value;
}

double* next_vector(FILE* json) {
  double* v = malloc(3*sizeof(double));
  expect_c(json, '[');
  skip_ws(json);
  v[0] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[1] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[2] = next_number(json);
  skip_ws(json);
  expect_c(json, ']');
  return v;
}

void read_scene(char* filename, Camera* camera, Object** objects) {
  int c;
  int current_item = -1; //for tracking the current object in the Object array
  int current_type; //for tracking the current object we are reading from the json list
  FILE* json = fopen(filename, "r");
  //if file does not exist
  if (json == NULL) {
    fprintf(stderr, "Error: Could not open file \"%s\"\n", filename);
    exit(1);
  }
  
  skip_ws(json);
  
  // Find the beginning of the list
  expect_c(json, '[');

  skip_ws(json);

  // Find the objects

  while (1) {
    c = fgetc(json);
    if (c == ']') {
      fprintf(stderr, "Error: This is the worst scene file EVER.\n");
      fclose(json);
      return;
    }
    if (c == '{') {
      skip_ws(json);
    
      // Parse the object
      char* key = next_string(json);
      if (strcmp(key, "type") != 0) {
        fprintf(stderr, "Error: Expected \"type\" key on line number %d.\n", line);
        exit(1);
      }

      skip_ws(json);

      expect_c(json, ':');

      skip_ws(json);

      char* value = next_string(json);

      if (strcmp(value, "camera") == 0) {
        current_type = 0;
      } 
      else if (strcmp(value, "sphere") == 0) {
        current_item++;
        if(current_item>127){
          fprintf(stderr, "Error: Too many objects in JSON. Program can only handle 128 objects.\n");
          exit(1);
        }
        objects[current_item] = malloc(sizeof(Object));
        objects[current_item]->type = 0;
        current_type = 1;
      } 
      else if (strcmp(value, "plane") == 0) {
        current_item++;
        if(current_item>127){
          fprintf(stderr, "Error: Too many objects in JSON. Program can only handle 128 objects.\n");
          exit(1);
        }
        objects[current_item] = malloc(sizeof(Object));
        objects[current_item]->type = 1;
        current_type = 2;
      } 
      else { 
        fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, line);
        exit(1);
      }
      skip_ws(json);

      while (1) {
        // , }
        c = next_c(json);
        if (c == '}') {
          // stop parsing this object
          break;
        } 
        else if (c == ',') {
          // read another field
          skip_ws(json);
          char* key = next_string(json);
          skip_ws(json);
          expect_c(json, ':');
          skip_ws(json);
          if (strcmp(key, "width") == 0){
            if(current_type == 0){  //only camera has width
              camera->width = next_number(json);
            }
            else{
              fprintf(stderr, "Error: Current object type has width value on line number %d.\n", line);
              exit(1);
            }
          }
          else if(strcmp(key, "height") == 0){
            if(current_type == 0){  //only camera has height
              camera->  height = next_number(json);
            }
            else{
              fprintf(stderr, "Error: Current object type has height value on line number %d.\n", line);
              exit(1);
            }
          }
          else if(strcmp(key, "radius") == 0){
            if(current_type == 1){  //only spheres have radius
              objects[current_item]->sphere.radius = next_number(json);
            }
            else{
              fprintf(stderr, "Error: Current object type cannot have radius value! Detected on line number %d.\n", line);
              exit(1);
            }
          }     
          else if(strcmp(key, "color") == 0){
            if(current_type == 1 || current_type == 2){  //only spheres and planes have color
                double* vector = next_vector(json);
                objects[current_item]->color[0] = (*vector);
                vector++;
                objects[current_item]->color[1] = (*vector);
                vector++;
                objects[current_item]->color[2] = (*vector); 
            }
            else{
              fprintf(stderr, "Error: Camera type has color value on line number %d.\n", line);
              exit(1);
            }
          } 
          else if(strcmp(key, "position") == 0){
            if(current_type == 1 || current_type == 2){  //only spheres and planes have position
              double* vector = next_vector(json);
              objects[current_item]->position[0] = (*vector++);
              objects[current_item]->position[1] = (*vector++);
              objects[current_item]->position[2] = (*vector++);  
            }
            else{
              fprintf(stderr, "Error: Camera type has position value on line number %d.\n", line);
              exit(1);
            }
          } 
          else if(strcmp(key, "normal") == 0){
            if(current_type == 2){  //only planes have normal
              double* vector = next_vector(json);
              objects[current_item]->plane.normal[0] = (*vector++);
              objects[current_item]->plane.normal[1] = (*vector++);
              objects[current_item]->plane.normal[2] = (*vector++);  
            }
            else{
              fprintf(stderr, "Error: Only planes have normal values on line number %d.\n", line);
              exit(1);
            }
          } 
          else{
            fprintf(stderr, "Error: Unknown property, \"%s\", on line %d.\n",
                key, line);
            //char* value = next_string(json);
          }

          skip_ws(json);
        } 
        else {
          fprintf(stderr, "Error: Unexpected value on line %d\n", line);
          exit(1);
        }
      }

      skip_ws(json);

      c = next_c(json);

      if (c == ',') {
        // noop
        skip_ws(json);
      } 
      else if (c == ']') {
        fclose(json);
        return;
      } 
      else {
        fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
        exit(1);
      }
    }
  }
}

double sqr(double v) {
  return v*v;
}

void normalize(double* v) {
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

double sphere_intersection(double* Ro, double* Rd, double* C, double r) {
  // Step 1. Find the equation for the object you are
  // interested in..  (e.g., sphere)
  //
  // x^2 + y^2 + z^2 = r^2
  //
  // Step 2. Parameterize the equation with a center point
  // if needed
  //
  // (x-Cx)^2 + (y-Cy)^2 + (z-Cz)^2 = r^2
  //
  // Step 3. Substitute the eq for a ray into our object
  // equation.
  //
  // (Rox + t*Rdx - Cx)^2 + (Roy + t*Rdy - Cy)^2 + (Roz + t*Rdz - Cz)^2 - r^2 = 0
  //
  // Step 4. Solve for t.
  //
  // Step 4a. Rewrite the equation (flatten).
  //
  // -r^2 +
  // t^2 * Rdx^2 +
  // t^2 * Rdy^2 +
  // t^2 * Rdz^2 +
  // 2*t * Rox * Rdx -
  // 2*t * Rdx * Cx +
  // 2*t * Roy * Rdy -
  // 2*t * Rdy * Cy +
  // 2*t * Roz * Rdz -
  // 2*t * Rdz * Cz +
  // Rox^2 -
  // 2*Rox*Cx +
  // Cx^2 +
  // Roy^2 -
  // 2*Roy*Cy +
  // Cy^2 +
  // Roz^2 -
  // 2*Roz*Cz +
  // Cz^2 = 0
  //
  // Steb 4b. Rewrite the equation in terms of t.
  //
  // t^2 * (Rdx^2 + Rdy^2 + Rdz^2) +
  // t * (2 * (Rox * Rdx - Rdx * Cx + Roy * Rdy - Rdy *Cy Roz * Rdz - Rdz * Cz)) +
  // Rox^2 - 2*Rox*Cx + Cx^2 + Roy^2 - 2*Roy*Cy + Cy^2  + Roz^2 - 2*Roz*Cz + Cz^2 - r^2 = 0
  //
  // Use the quadratic equation to solve for t..
  double a = (sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]));
  double b = (2 * (Ro[0] * Rd[0] - Rd[0] * C[0] + Ro[1] * Rd[1] - Rd[1] * C[1] + Ro[2] * Rd[2] - Rd[2] * C[2]));
  double c = sqr(Ro[0]) - 2*Ro[0]*C[0] + sqr(C[0]) + sqr(Ro[1]) - 2*Ro[1]*C[1] + sqr(C[1]) + sqr(Ro[2]) - 2*Ro[2]*C[2] + sqr(C[2]) - sqr(r);

  double det = sqr(b) - 4 * a * c;
  if (det < 0) return -1;

  det = sqrt(det);

  double t0 = (-b - det) / (2*a);
  if (t0 > 0) return t0;

  double t1 = (-b + det) / (2*a);
  if (t1 > 0) return t1;

  return -1;
}

double plane_intersection(double* Ro, double* Rd, double* P, double* N) {
  // Step 1. Find the equation for the object you are
  // interested in..  (e.g., plNW)
  //
  // Ax + By + Cz + D = 0 
  //
  // Step 2. Parameterize the equation with a center point
  // if needed
  //
  // Since N = [A, B, C]
  //
  // Nx(x-Px) + Ny(y-Py) + Nz(z-Pz) = 0
  //
  // Step 3. Substitute the eq for a ray into our object
  // equation.
  //
  // Nx(Rox + t*Rdx - Px) + Ny(Roy + t*Rdy - Py) + Nz(Roz + t*Rdz - Pz) = 0
  //
  // Step 4. Solve for t.
  //
  // Step 4a. Rewrite the equation (flatten).
  //
  // NxRox + t*Nx*Rdx - NxPx + NyRoy + t*Ny*Rdy - NyPy + NzRoz + t*Nz*Rdz - NzPz = 0
  // 
  // 
  //
  // Steb 4b. Rewrite the equation in terms of t.
  //
  // t(Nx*Rdx + Ny*Rdy + Nz*Rdz) = Nx*Px + Ny*Py + Nz*Pz - Nx*Rox - Ny*Roy - Nz*Roz
  //
  // t = (NxPx + NyPy + NzPz - NxRox - NyRoy - NzRoz)/(Nx*Rdx + Ny*Rdy + Nz*Rdz) 
  //
  //normalize the normal vector
  normalize(N);
  double t = (N[0]*P[0] + N[1]*P[1] + N[2]*P[2] - N[0]*Ro[0] - N[1]*Ro[1] - N[2]*Ro[2])/(N[0]*Rd[0] + N[1]*Rd[1] + N[2]*Rd[2]); 
  if (t > 0) return t;

  return -1;
}

void generate_scene(Camera* camera, Object** objects, Pixel* buffer, int width, int height){
  //write these objects to a ppm image
  double camera_width = camera->width;
  double camera_height = camera->height;
  double pixheight = camera_height / height;
  double pixwidth = camera_width / width;
  for (int y = 0; y < height; y += 1) {
    for (int x = 0; x < width; x += 1) {
      double Ro[3] = {0, 0, 0};
      // Rd = normalize(P - Ro)
      double Rd[3] = {
        0 - (camera_width/2) + pixwidth * (x + 0.5),
        0 - (camera_height/2) + pixheight * (y + 0.5),
        1
      };
      normalize(Rd);

      double best_t = INFINITY;
      for (int i=0; objects[i] != 0; i += 1) {
        double t = 0;

        switch(objects[i]->type) {
        case 0:
          t = sphere_intersection(Ro, Rd, objects[i]->position, objects[i]->sphere.radius);
          break;
        case 1:
          t = plane_intersection(Ro, Rd, objects[i]->position, objects[i]->plane.normal);
          break;
        default:
          // Horrible error
          exit(1);
        }
        if (t > 0 && t < best_t) {
          best_t = t;
          int position = (height-(y+1))*width+x;
          buffer[position].r = (unsigned char) 255*objects[i]->color[0];
          buffer[position].g = (unsigned char) 255*objects[i]->color[1];
          buffer[position].b = (unsigned char) 255*objects[i]->color[2];
        }
        if (best_t > 0 && best_t != INFINITY) {

        } 
        else {
          int position = (height-(y+1))*width+x;
          buffer[position].r = 0;
          buffer[position].g = 0;
          buffer[position].b = 0;
        }
      }
    }
  } 
}

void write_p3(Pixel *buffer, FILE *output_file, int width, int height, int max_color){
  fprintf(output_file, "P3\n%d %d\n%d\n", width, height, max_color);
  int current_width = 1;
  for(int i = 0; i < width*height; i++){
    fprintf(output_file, "%d %d %d ", buffer[i].r, buffer[i].g, buffer[i].b);
    if(current_width >= 70%12){ //ppm line length = 70, max characters to pixels = 12.
      fprintf(output_file, "\n");
      current_width = 1;
    }
    else{
      current_width++;
    }
  }
}