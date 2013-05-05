(declare (usual-integrations))

(define (flo:sigma f low high)
  (if (fix:> low high)
    0.
    (let lp ((i (fix:+ low 1)) (sum (f low)))
      (if (fix:> i high)
        sum
        (lp (fix:+ i 1) (flo:+ sum (f i)))))))

(define (flo:make-initialized-vector n proc)
  (let ((result (flo:vector-cons n)))
    (let lp ((i 0))
      (if (fix:< i n)
        (begin
          (flo:vector-set! result i (proc i))
          (lp (fix:+ i 1)))
        result))))

(define (list->flo-vector lst)
  (let* ((len (length lst))
         (v (flo:vector-cons len)))
    (let lp ((lst lst)
             (idx 0))
      (if (< idx len)
        (begin
          (flo:vector-set! v idx (exact->inexact (car lst)))
          (lp (cdr lst) (+ idx 1)))
        v))))

(define (vector->flo-vector vec)
  (let* ((len (vector-length vec))
         (v (flo:vector-cons len)))
    (let lp ((idx 0))
      (if (< idx len)
        (begin
          (flo:vector-set! v idx (exact->inexact (vector-ref vec idx)))
          (lp (+ idx 1)))
        v))))

(define (pair->flo-vector pair)
  (let ((v (flo:vector-cons 2)))
    (flo:vector-set! v 0 (exact->inexact (car pair)))
    (flo:vector-set! v 1 (exact->inexact (cdr pair)))
    v))

(define (flo:vector-sum v)
  (let ((len (flo:vector-length v)))
    (let lp ((idx 0)
             (tot 0.))
      (if (< idx len)
        (lp (+ idx 1) (flo:+ tot (flo:vector-ref v idx)))
        tot))))

(define (flo:sum . args)
  (let lp ((tot 0.)
           (lst args))
    (if (null? lst)
      tot
      (lp (flo:+ tot (car lst)) (cdr lst)))))

