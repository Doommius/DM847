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

#ifndef PICKING_H
#define PICKING_H

#include "ims_core.h"

#define picking(func_name) IMSPeakList* func_name(IMSPeakList* peak_list, pmap* pipeline_parameters)

typedef IMSPeakList* (*picking_p)(IMSPeakList*, pmap* pipeline_parameters);

picking(merging_by_signal);
picking(cluster_editing);
picking(em_clustering);
picking(no_picking);

double sim_score(IMSPeak* pointOne, IMSPeak* pointTwo, pmap* pipeline_parameters);

#endif /* PICKING_H */
