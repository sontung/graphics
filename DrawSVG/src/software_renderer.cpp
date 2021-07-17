#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {

    // set top level transformation
    transformation = svg_2_screen;

    // draw all elements
    for ( size_t i = 0; i < svg.elements.size(); ++i ) {
        draw_element(svg.elements[i]);
    }

    // draw canvas outline
    Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
    Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
    Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
    Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

    rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
    rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
    rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
    rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

    // resolve and send to render target
    resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

    // Task 4:
    // You may want to modify this for supersampling support
    this->sample_rate = sample_rate;
    this->sample_buffer_w = this->target_w*sample_rate;
    this->sample_buffer_h = this->target_h*sample_rate;
    sample_buffer.resize(4 * sample_buffer_w * sample_buffer_h);
    counting_buffer.resize(4*target_w*target_h);
    memset(&sample_buffer[0], 0.0, 4 * sample_buffer_w * sample_buffer_h*sizeof(float));
    memset(&counting_buffer[0], 0, target_w * target_h*sizeof(int));

}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

    // Task 4:
    // You may want to modify this for supersampling support
    this->render_target = render_target;
    this->target_w = width;
    this->target_h = height;

    this->sample_buffer_w = width*sample_rate;
    this->sample_buffer_h = height*sample_rate;
    sample_buffer.resize(4 * sample_buffer_w * sample_buffer_h);
    counting_buffer.resize(4*target_w*target_h);
    memset(&sample_buffer[0], 0.0, 4 * sample_buffer_w * sample_buffer_h*sizeof(float));
    memset(&counting_buffer[0], 0, target_w * target_h*sizeof(int));

}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

    // Task 5 (part 1):
    // Modify this to implement the transformation stack

    switch(element->type) {
    case POINT:
        draw_point(static_cast<Point&>(*element));
        break;
    case LINE:
        draw_line(static_cast<Line&>(*element));
        break;
    case POLYLINE:
        draw_polyline(static_cast<Polyline&>(*element));
        break;
    case RECT:
        draw_rect(static_cast<Rect&>(*element));
        break;
    case POLYGON:
        draw_polygon(static_cast<Polygon&>(*element));
        break;
    case ELLIPSE:
        draw_ellipse(static_cast<Ellipse&>(*element));
        break;
    case IMAGE:
        draw_image(static_cast<Image&>(*element));
        break;
    case GROUP:
        draw_group(static_cast<Group&>(*element));
        break;
    default:
        break;
    }

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

    Vector2D p = transform(point.position);
    rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

    Vector2D p0 = transform(line.from);
    Vector2D p1 = transform(line.to);
    rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

    Color c = polyline.style.strokeColor;

    if( c.a != 0 ) {
        int nPoints = polyline.points.size();
        for( int i = 0; i < nPoints - 1; i++ ) {
            Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
            Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
            rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
        }
    }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

    Color c;

    // draw as two triangles
    float x = rect.position.x;
    float y = rect.position.y;
    float w = rect.dimension.x;
    float h = rect.dimension.y;

    Vector2D p0 = transform(Vector2D(   x   ,   y   ));
    Vector2D p1 = transform(Vector2D( x + w ,   y   ));
    Vector2D p2 = transform(Vector2D(   x   , y + h ));
    Vector2D p3 = transform(Vector2D( x + w , y + h ));

    // draw fill
    c = rect.style.fillColor;
    if (c.a != 0 ) {
        rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
        rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
    }

    // draw outline
    c = rect.style.strokeColor;
    if( c.a != 0 ) {
        rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
        rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
        rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
        rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
    }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

    Color c;

    // draw fill
    c = polygon.style.fillColor;
    if( c.a != 0 ) {

        // triangulate
        vector<Vector2D> triangles;
        triangulate( polygon, triangles );

        // draw as triangles
        for (size_t i = 0; i < triangles.size(); i += 3) {
            Vector2D p0 = transform(triangles[i + 0]);
            Vector2D p1 = transform(triangles[i + 1]);
            Vector2D p2 = transform(triangles[i + 2]);
            rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
        }
    }

    // draw outline
    c = polygon.style.strokeColor;
    if( c.a != 0 ) {
        int nPoints = polygon.points.size();
        for( int i = 0; i < nPoints; i++ ) {
            Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
            Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
            rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
        }
    }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

    // Extra credit

}

void SoftwareRendererImp::draw_image( Image& image ) {

    Vector2D p0 = transform(image.position);
    Vector2D p1 = transform(image.position + image.dimension);

    rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

    for ( size_t i = 0; i < group.elements.size(); ++i ) {
        draw_element(group.elements[i]);
    }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

    // fill in the nearest pixel
    int sx = (int) floor(x);
    int sy = (int) floor(y);

    // check bounds
    if ( sx < 0 || sx >= target_w ) return;
    if ( sy < 0 || sy >= target_h ) return;

    // fill sample - NOT doing alpha blending!
    render_target[4 * (sx + sy * target_w)    ] = (uint8_t) (color.r * 255);
    render_target[4 * (sx + sy * target_w) + 1] = (uint8_t) (color.g * 255);
    render_target[4 * (sx + sy * target_w) + 2] = (uint8_t) (color.b * 255);
    render_target[4 * (sx + sy * target_w) + 3] = (uint8_t) (color.a * 255);

}

void swap(float &u, float &v) {
    float temp = u;
    u=v;
    v=temp;
}

float fpart(float x) {
    return x-floor(x);
}

float rfpart(float x) {
    return 1-fpart(x);
}

float round(float x) {
    return floor(x+0.5);
}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

    // Task 2:
    // Implement line rasterization
    bool steep = abs(y1-y0)>abs(x1-x0);
    if (steep) {
        swap(x0, y0);
        swap(x1, y1);
    }
    if (x0 > x1) {
        swap(x0, x1);
        swap(y0, y1);
    }
    float dx = x1-x0;
    float dy = y1-y0;
    float gradient;
    if (abs(dx)<0.0001) gradient = 1.0;
    else gradient=dy/dx;

    float xend=round(x0);
    float yend=y0+gradient*(xend-x0);
    float xgap=rfpart(x0+0.5);
    float xpxl1 = xend;
    float ypxl1 = floor(yend);
    if (steep) {
        rasterize_point(ypxl1,   xpxl1, rfpart(yend) * xgap*color);
        rasterize_point(ypxl1+1, xpxl1,  fpart(yend) * xgap*color);
    } else {
        rasterize_point(xpxl1, ypxl1  , rfpart(yend) * xgap*color);
        rasterize_point(xpxl1, ypxl1+1,  fpart(yend) * xgap*color);
    }
    float intery = yend+gradient;

    xend=round(x1);
    yend=y1+gradient*(xend-x1);
    xgap=fpart(x1+0.5);
    float xpxl2=xend;
    float ypxl2=floor(yend);
    if (steep) {
        rasterize_point(ypxl2  , xpxl2, rfpart(yend) * xgap*color);
        rasterize_point(ypxl2+1, xpxl2,  fpart(yend) * xgap*color);
    } else {
        rasterize_point(xpxl2, ypxl2,  rfpart(yend) * xgap*color);
        rasterize_point(xpxl2, ypxl2+1, fpart(yend) * xgap*color);
    }
    if (steep) {
        for (float x=xpxl1+1; x<=xpxl2-1; x++) {
            rasterize_point(floor(intery), x, rfpart(intery)*color);
            rasterize_point(floor(intery)+1, x, fpart(intery)*color);
            intery = intery+gradient;
        }
    } else {
        for (float x=xpxl1+1; x<=xpxl2-1; x++) {
            rasterize_point(x, floor(intery), rfpart(intery)*color);
            rasterize_point(x, floor(intery)+1, fpart(intery)*color);
            intery = intery+gradient;
        }
    }
}

float point_line_realtion(float px, float py,
                            float x0, float y0,
                            float x1, float y1){
  float dx = x1- x0;
  float dy = y1 - y0;
  float dpx = px - x0;
  float dpy = py - y0;
  float relation = -dpx * dy + dpy * dx;
  return relation;
}

bool point_in_triangle(float px, float py,
                      float x0, float y0,
                      float x1, float y1,
                      float x2, float y2){
  float r1 = point_line_realtion(px, py, x0, y0, x1, y1);
  float r2 = point_line_realtion(px, py, x1, y1, x2, y2);
  float r3 = point_line_realtion(px, py, x2, y2, x0, y0);
  if ((r1*r2 < 0) || (r2 * r3 < 0)) return false;
  return true;
}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
    // Task 3:
    // Implement triangle rasterization
    float sample_step = 1.0/ static_cast<float>(sample_rate);

    float x_start = floor(min({x0, x1, x2})-1)+sample_step/2;
    float y_start_mark = floor(min({y0, y1, y2})-1)+sample_step/2;
    float x_end = floor(max({x0, x1, x2})+1)+sample_step/2;
    float y_end = floor(max({y0, y1, y2})+1)+sample_step/2;
    float y_start;
    if (x_start<sample_step/2) x_start=sample_step/2;
    if (x_end>target_w) x_end=target_w;
    if (y_start_mark<sample_step/2) y_start_mark=sample_step/2;
    if (y_end>target_h) y_end=target_h;
    while (x_start<x_end) {
        y_start=y_start_mark;
        while (y_start<y_end) {
            if (point_in_triangle(x_start, y_start, x0, y0, x1, y1, x2, y2)) {
                fill_sample(x_start, y_start, color);
                rasterize_point(x_start, y_start, color);
            }
            y_start+=sample_step;
        }
        x_start += sample_step;
    }
}

void SoftwareRendererImp::fill_sample( float x, float y, const Color& color ) {
    // fill in the nearest pixel
    int sx;
    int sy;

    // fill countinng buffer
    sx = (int) floor(x);
    sy = (int) floor(y);
    if ( sx < 0 || sx >= target_w ) return;
    if ( sy < 0 || sy >= target_h ) return;
    counting_buffer[sx + sy * target_w]++;

    return;
}

void SoftwareRendererImp::fill_pixel(int x, int y) {

    int r_old = static_cast<int>(render_target[4 * (x + y * target_w)    ]);
    int g_old = static_cast<int>(render_target[4 * (x + y * target_w) + 1]);
    int b_old = static_cast<int>(render_target[4 * (x + y * target_w) + 2]);
    int a_old = static_cast<int>(render_target[4 * (x + y * target_w) + 3]);
    int weight = counting_buffer[x + y * target_w];

    if (weight > 0) {
        float full = static_cast<float>(sample_rate*sample_rate);
        render_target[4 * (x + y * target_w)    ] = (uint8_t) (r_old*weight/full);
        render_target[4 * (x + y * target_w) + 1] = (uint8_t) (g_old*weight/full);
        render_target[4 * (x + y * target_w) + 2] = (uint8_t) (b_old*weight/full);
        render_target[4 * (x + y * target_w) + 3] = (uint8_t) (a_old*weight/full);
    }
    return;
}


void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
    // Task 6:
    // Implement image rasterization

}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

    // Task 4:
    // Implement supersampling
    // You may also need to modify other functions marked with "Task 4".
    cout<<"resolving \n"<<endl;
    for (int x=0; x<target_w; x++) for (int y=0; y<target_h; y++) fill_pixel(x, y);
    memset(&counting_buffer[0], 0, target_w * target_h*sizeof(int));
    return;

}


} // namespace CMU462
