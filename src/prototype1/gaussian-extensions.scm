;; TODO remove these
(define default:+ +)
(define defualt:* *)

;; apply this method if there are unforced gaussians in args
(define (gaussian:+ . args)
  ((partition gaussian-sample? args)
   (lambda (gaussians nongaussians)
     (gaussian-builder:new-sum gaussians (apply default:+ (map force-sample nongaussians))))))

;; apply this method if there's one unforced sample and it's a gaussian
(define (gaussian:* . args)
  ((partition gaussian-sample? args)
   (lambda (gaussians nongaussians)
     (if (not (= (length gaussians) 1))
       (error "BRO"))
     (gaussian-builder:new-scaling (car gaussians) (apply default:* (map force-sample nongaussians))))))



;; new gaussian object, not a generic, hooks into the joint
(define (gaussian:new mean var)
  (let ((var (force-sample var)))
    (if (gaussian-sample? mean)
      (gaussian:+ mean (gaussian:new 0 var))
      (gaussian-builder:new-independent mean var))))


;; Gaussian Builder ;;

(define *cov-factor* (make-vector 16))
(define *sources-list* '())
(define *cov-factor-size* 0)
(define *indices* (make-eq-hash-table))

(load "randutil")
(load "lazy-collapsed-gaussians")

;; (gaussian-builder:new-independent mean var)
;;   makes and returns a new gausisan sample object corresponding within the
;;   joint gausisan to an independent sample (a "source")
(define (gaussian-builder:new-independent mean var)
  (let ((sample-obj (%gaussian-sample:new (gaussian:make-params mean var) 'the-unsampled-value)))
    (let ((index (gaussian-builder:get-index sample-obj)))
      (set! *sources-list* (cons sample-obj *sources-list*))
      (gaussian-builder:add-row! `#(,(cons index (sqrt var))))
      sample-obj)))

;; (gaussian-builder:new-sum list-of-unsorted-existing-gaussians number)
;;   makes and returns a new gausisan sample object correspoding within the
;;   joint gaussian to a sum of the given gaussians and the numerical constant
(define (gaussian-builder:new-sum list-of-gaussians const)
  (define (sum-squares sparse-row-vec)
    (apply default:+ (map square (vector->list (vector-map cdr sparse-row-vec)))))
  (define (sum-of-means gaussian-objects)
    (apply default+ (map gaussian:mean (map gaussian-sample:params gaussian-objects))))

  (let ((newrow (gaussian-builder:row-dot list-of-gaussians)))
    (let ((newvar (sum-squares newrow))
          (newmean (sum-of-means list-of-gaussians)))
      (let* ((sample-obj (%gaussian-sample:new (gaussian:make-params newmean newvar) 'the-unsampled-value))
             (index (gaussian-builder:get-index sample-obj)))
        (vector-set! newrow 0 (cons index newvar))
        (gaussian-builder:add-row! newrow)
        sample-obj))))

;; (gaussian-builder:new-scaling existing-gaussian-obj number)
;;   makes and returns a new gaussian sample object corresponding within the
;;   joint gaussian to a scaling of the given gaussian by the numerical constant
;; (define (gaussian-builder:new-scaling gaussian const)
;;   )


;; (define (gaussian-builder:get-final-params)
;;   )

;; (define (gaussian-builder:sample!)
;;   )

;; instead of all this i could just do dense matrix ops to begin with

(define (gaussian-builder:add-row! row)
  (let ((len (vector-length *cov-factor*)))
    (if (= *cov-factor-size* len)
      (set! *cov-factor* (vector-grow *cov-factor* (* 2 len)))))
  (vector-set! *cov-factor* *cov-factor-size* row)
  (set! *cov-factor-size* (default:+ *cov-factor-size* 1)))

(define (gaussian-builder:get-index obj)
  (hash-table/lookup *indices* obj
                     (lambda (idx) idx)
                     (lambda ()
                       (let ((new-index (hash-table/count *indices*)))
                         (hash-table/put! *indices* obj new-index)
                         new-index))))

(define (gaussian-builder:row-dot list-of-gaussians)
  (let ((result-hash (make-eq-hash-table)))
    (let lp1 ((lst list-of-gaussians))
      (if (not (null? lst))
        (let* ((row (vector-ref *cov-factor* (gaussian-builder:get-index (car lst))))
               (len (vector-length row)))
          (let lp2 ((idx 0))
            (if (< idx len)
              (let* ((pair (vector-ref row idx))
                     (col (car pair))
                     (val (cdr pair)))
                (hash-table/lookup result-hash
                                   col
                                   (lambda (oldval)
                                     (hash-table/put! result-hash col (default:+ val oldval)))
                                   (lambda ()
                                     (hash-table/put! result-hash col val)))
                (lp2 (default:+ idx 1)))))
          (lp1 (cdr lst)))))
    (list->vector (cons (cons #f 0) (hash-table->alist result-hash)))))


;; NOTE: choice vals could be samples, in which case set prior score will fail.
;; prior score should be computed on the fly instead of being stored! TODO

;; one big emit should definitely be supported
;; what about multiple emits? just do the environment thing, require every
;; execution exactly explains all the emitted data, then this system should
;; work on its own (with emit just setting values and the final emit, which
;; knows it's final because it explains the last datapoint, and then doing
;; usual thing)

;; emit only sets values that are uncollapsed; conditioning on a collapsed
;; value means it needs a likelihood, which could be provided by an uncollapsed

;; NOTE: could have a treewidth measure (largest number of g terms in a sum)
;; and do message passing below some threshold
