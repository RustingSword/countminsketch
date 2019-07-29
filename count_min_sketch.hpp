#ifndef COUNT_MIN_SKETCH_H
#define COUNT_MIN_SKETCH_H
/**
  Daniel Alabi
  Count-Min Sketch Implementation based on paper by
  Muthukrishnan and Cormode, 2004
 **/

#include <inttypes.h>
#include <vector>
// define some constants
# define LONG_PRIME 4294967311ul

/** CountMinSketch class definition here **/
class CountMinSketch {
  // width, depth
  uint32_t w,d;

  // eps (for error), 0.01 < eps < 1
  // the smaller the better
  float eps;

  // gamma (probability for accuracy), 0 < gamma < 1
  // the bigger the better
  float gamma;

  // aj, bj \in Z_p
  // both elements of fild Z_p used in generation of hash
  // function
  unsigned int aj, bj;

  // total count so far
  uint32_t total;

  // array of arrays of counters
  uint32_t **C;

  // array of hash values for a particular item
  // contains two element arrays {aj,bj}
  uint64_t **hashes;

  // generate "new" aj,bj
  void genajbj(uint64_t **hashes, int i);

  public:
  // constructor
  CountMinSketch(float eps = 0.01, float gamma = 0.01);

  // update item (int) by count c
  void update(uint64_t item, int c = 1);
  // update item (string) by count c
  void update(const char *item, int c = 1);

  // estimate count of item i and return count
  uint32_t estimate(uint64_t item);
  uint32_t estimate(const char *item);

  // return total count
  uint32_t totalcount();

  // generates a hash value for a string
  // same as djb2 hash function
  uint64_t hashstr(const char *str);

  // clear all counters
  void reset();

  // destructor
  ~CountMinSketch();
};

/*
 * \brief circular CMS, can be used on streaming data to record item count of a
 * specified time span, e.g. 1 hour.
 *
 * currently for integers only.
 */
class CircularCMS {
  public:
    explicit CircularCMS(int num_counter) :
      num_counter_(num_counter),
      current_counter_index_(-1) {
        for (int i = 0; i < num_counter_; ++i) {
          counters_.push_back(new CountMinSketch());
        }
      }
    ~CircularCMS() {
      for (size_t i = 0; i < counters_.size(); ++i)
        delete counters_[i];
    }
    int update(uint64_t item, int counter_index, int c = 1) {
      if (counter_index > static_cast<int>(counters_.size())) {
        std::cerr << "[ERROR] counter index overflow (" << counter_index
                  << " >= " << counters_.size() << std::endl;
        return -1;
      }
      if (counter_index != current_counter_index_) {
        current_counter_index_ = counter_index;
        counters_[current_counter_index_]->reset();
      }
      counters_[current_counter_index_]->update(item, c);
      return 0;
    }
    int estimate(uint64_t item) {
      int count = 0;
      for (size_t i = 0; i < counters_.size(); ++i)
        count += counters_[i]->estimate(item);
      return count;
    }

  private:
    std::vector<CountMinSketch*> counters_;
    int num_counter_;
    int current_counter_index_;
};

/* \brief Use multiple CircularCMS together in multi-thread scenario
*/
class MultiCMS {
  public:
    explicit MultiCMS(int num_thread, int count_interval_seconds) :
      num_thread_(num_thread), count_interval_seconds_(count_interval_seconds) {
        if (num_thread < 0 || count_interval_seconds < 0 ||
            count_interval_seconds > 3600) {
          std::cerr << "[ERROR] wrong num_thread (" << num_thread
                    << ") or count_interval_seconds (" << count_interval_seconds
                    << ")\n";
          return;
        }
        num_counter_ = int(ceil(3600.0 / count_interval_seconds));
        for (int i = 0; i < num_thread; ++i) {
          counters_.push_back(new CircularCMS(num_counter_));
        }
      }
    ~MultiCMS() {
      for (size_t i = 0; i < counters_.size(); ++i) {
        delete counters_[i];
      }
    }
    int update(int thread_id, uint64_t item, int c = 1) {
      time_t ts = time(NULL);
      int counter_index = (ts % 3600) / count_interval_seconds_;
      if ((thread_id < 0) || (thread_id >= num_thread_) ||
          (counter_index >= num_counter_)) {
        std::cerr << "[ERROR] wrong thread id (" << thread_id
                  << ") or counter_index (" << counter_index << ")\n";
        return -1;
      }
      return counters_[thread_id]->update(item, counter_index, c);
    }
    int estimate(uint64_t item) {
      int count = 0;
      for (size_t i = 0; i < counters_.size(); ++i) {
        count += counters_[i]->estimate(item);
      }
      return count;
    }

  private:
    int num_counter_;
    int num_thread_;
    int count_interval_seconds_;
    std::vector<CircularCMS*> counters_;
};

#endif  // COUNT_MIN_SKETCH_H
