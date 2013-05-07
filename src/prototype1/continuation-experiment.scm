
(define global 0)

(define (>>= ma a->mb ) 
  (lambda (in-state)
    (let ((a-&-state (ma in-state))
          (mb (f (car a))))
      (mb (cdr a-&-state)))))

(define (>> ma mb)
  (>>= ma (lambda (x) mb)))

(define (return x)
  (lambda (s)
    (pair x s)))

(define (sample global)
  (call-with-current-continuation
    (lambda (k)
      (fluid-let ((global (+ global 1)))
        (pp global)
        (k (cons 7 global))))))

(define *continuations* '())
(define (get-previous-continuation index)
  (if (> index (- (length *continuations*) 1))
      (error "Index greater than number of stored continuations in *continuations*")
      (let ((j (- (length *continuations*) index 1)))
        (set! *continuations* (sublist *continuations* j (length *continuations*)))
        (car *continuations*))))

(define (sample2)
  (let ((current-global global))
    (call-with-current-continuation
     (lambda (k)
       (let ((new-k (lambda (val)
                     (begin 
                       (set! global current-global)
                       (k val)))))
         (set! *continuations* (cons new-k *continuations*))
         (set! global (+ 1 current-global))
         (pp global)
         (k 7))))))
       
(define (user-code)
  (let ((x (sample2)))
    (let ((y (sample2)))
      (let ((z (sample2)))
      z))))


  
