/*!
 * @author Marianna D'Addario, Dominik Kopczynski
 * @e-mail {Marianna.Daddario, Dominik.Kopczynski}@tu-dortmund.de
 * Copyright (c) 2013 Marianna D'Addario, Dominik Kopczynski
 *
 * LICENSE:
 * PEAX is a free software for academic use only.
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef CANDIDATE_DETECTION_H
#define CANDIDATE_DETECTION_H

#include "ims_core.h"

#define X_ALIGN 0
#define Y_ALIGN 1
#define GAP 0.1

#define DEF_AREASIZE 10
#define DEF_MINPEAKSIGNAL 5
#define DEF_MARGINSIZE 2

#define candidate_detection(func_name) IMSPeakList* func_name(IMSMeasurement* measurement, pmap* pipeline_parameters)

typedef IMSPeakList* (*candidate_detection_p)(IMSMeasurement*, pmap* pipeline_parameters);

struct point{
    int x;
    int y;
};

struct line{
    point ll; // lower left
    point ur; // upper right
    vector< point* >* points;
};


struct align{
    int old_index;
    int new_index;
};

candidate_detection(local_maxima);
candidate_detection(cross_finding);
candidate_detection(bit_parallel);


#endif /* CANDIDATE_DETECTION_H */
